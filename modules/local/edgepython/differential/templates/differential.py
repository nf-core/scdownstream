#!/usr/bin/env python3

import platform
import re
from pathlib import Path

import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd
import patsy
import yaml

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgePython pseudobulk differential expression")

STANDARD_DE_COLUMNS = ["gene", "log2fc", "pvalue", "padj", "group", "contrast", "stratum"]


def write_standard_de_parquet(
    df,
    path,
    *,
    gene_col,
    log2fc_col,
    pvalue_col,
    padj_col,
    group,
    contrast,
    stratum="",
):
    out = df.copy()
    if gene_col is None:
        if "genes" in out.columns:
            gene_col = "genes"
        else:
            out = out.reset_index()
            gene_col = out.columns[0]

    rename = {gene_col: "gene", log2fc_col: "log2fc"}
    if pvalue_col and pvalue_col in out.columns:
        rename[pvalue_col] = "pvalue"
    if padj_col and padj_col in out.columns:
        rename[padj_col] = "padj"
    out = out.rename(columns=rename)

    if "pvalue" not in out.columns:
        out["pvalue"] = pd.NA
    if "padj" not in out.columns:
        out["padj"] = pd.NA

    n_rows = len(out)
    out["gene"] = out["gene"].astype(str)
    out["log2fc"] = pd.to_numeric(out["log2fc"], errors="coerce")
    out["pvalue"] = pd.to_numeric(out["pvalue"], errors="coerce")
    out["padj"] = pd.to_numeric(out["padj"], errors="coerce")

    def _meta_column(value):
        if isinstance(value, pd.Series):
            values = value.astype(str).to_numpy()
            if len(values) != n_rows:
                raise ValueError("Metadata length does not match table rows")
            return values
        return "" if value is None else str(value)

    out["group"] = _meta_column(group)
    out["contrast"] = _meta_column(contrast)
    out["stratum"] = _meta_column(stratum)

    extra = [column for column in out.columns if column not in STANDARD_DE_COLUMNS]
    out = out[STANDARD_DE_COLUMNS + extra]
    # Absolute path so pyarrow does not treat colons in Nextflow prefixes as URI schemes
    out.to_parquet(str(Path(path).resolve()), index=False, engine="pyarrow")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))


def donors_span_conditions(metadata: pd.DataFrame) -> bool:
    """True when at least one donor appears in more than one condition (paired design)."""
    conditions_per_donor = metadata.groupby("donor", observed=True)["condition"].apply(
        lambda values: values.astype(str).nunique()
    )
    return bool((conditions_per_donor > 1).any())


def condition_coef_column(design: pd.DataFrame, treatment: str) -> str | None:
    for name in (
        f"condition{treatment}",
        f"condition[T.{treatment}]",
        f"condition[{treatment}]",
    ):
        if name in design.columns:
            return name
    for name in design.columns:
        if name.startswith("condition") and treatment in name:
            return name
    return None


required_cols = ["donor", "condition", "celltype"]
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in pseudobulk adata.obs")

written = []
for celltype, celltype_data in adata.obs.groupby("celltype", observed=True):
    sub = adata[celltype_data.index].copy()
    metadata = sub.obs[["donor", "condition"]].copy()
    conditions = sorted(metadata["condition"].astype(str).unique())
    donors = sorted(metadata["donor"].astype(str).unique())
    if len(conditions) < 2 or len(donors) < 2:
        continue

    treatments = [condition for condition in conditions if condition != reference_condition]
    if not treatments:
        continue

    counts = sub.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    counts = np.asarray(counts).T

    metadata = metadata.copy()
    metadata["condition"] = metadata["condition"].astype(str)
    metadata["donor"] = metadata["donor"].astype(str)
    metadata["condition"] = pd.Categorical(
        metadata["condition"],
        categories=[reference_condition] + treatments,
    )

    design_formula = "~ donor + condition" if donors_span_conditions(metadata) else "~ condition"
    # edgepython.model_matrix drops coefficient names; keep a patsy DataFrame for lookup.
    design_df = patsy.dmatrix(design_formula, data=metadata, return_type="dataframe")
    design = np.asarray(design_df, dtype=np.float64)
    y = ep.make_dgelist(
        counts=counts,
        samples=metadata,
        genes=pd.DataFrame({"genes": np.asarray(sub.var_names.astype(str))}),
    )
    y = ep.calc_norm_factors(y)
    y = ep.estimate_disp(y)
    fit = ep.glm_ql_fit(y, design)
    safe_celltype = safe_name(celltype)

    for treatment in treatments:
        coef_name = condition_coef_column(design_df, treatment)
        if coef_name is None:
            continue
        coef_index = list(design_df.columns).index(coef_name)
        res = ep.glm_ql_ftest(fit, coef=coef_index)
        top = ep.top_tags(res, n=sub.n_vars)
        out_path = f"{prefix}_{safe_celltype}_{safe_name(treatment)}_results.parquet"
        write_standard_de_parquet(
            top["table"],
            out_path,
            gene_col=None,
            log2fc_col="logFC",
            pvalue_col="PValue",
            padj_col="FDR",
            group=treatment,
            contrast=f"{treatment} vs {reference_condition}",
            stratum=f"celltype={celltype}",
        )
        written.append(out_path)

if not written:
    raise ValueError("No cell-type strata passed replicate filtering for edgePython pseudobulk")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
