#!/usr/bin/env python3

import os
import platform
import re
from pathlib import Path

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
donor_col = "${donor_col}"
condition_col = "${condition_col}"
celltype_col = "${celltype_col}"
celltype_value = "${celltype_value}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgepython_sc differential expression")

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


for col in [donor_col, condition_col, celltype_col]:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")

subset = adata[adata.obs[celltype_col].astype(str) == str(celltype_value)].copy()
if subset.n_obs < 10:
    raise ValueError(f"Too few cells ({subset.n_obs}) in cell type '{celltype_value}'")

condition_values = subset.obs[condition_col].astype(str)
conditions = sorted(condition_values.unique())
donors = sorted(subset.obs[donor_col].astype(str).unique())
if len(conditions) < 2:
    raise ValueError("At least two conditions are required for differential expression")
if len(donors) < 2:
    raise ValueError("At least two donors are required for edgepython_sc")

treatments = [condition for condition in conditions if condition != reference_condition]
if not treatments:
    raise ValueError(f"Reference condition '{reference_condition}' is the only condition present")

counts = subset.X
if sp.issparse(counts):
    counts = counts.toarray()
counts = np.asarray(counts).T

design = pd.DataFrame(
    {"Intercept": np.ones(subset.n_obs, dtype=float)},
    index=subset.obs_names,
)
for treatment in treatments:
    design[treatment] = (condition_values == treatment).astype(float)

sample_ids = subset.obs[donor_col].astype(str).to_numpy()

fit = ep.glm_sc_fit(
    counts,
    design=design,
    sample=sample_ids,
    norm_method="TMM",
)
fit = ep.shrink_sc_disp(fit, robust=True)

written = []
gene_mask = fit.get("gene_mask")
if gene_mask is not None:
    gene_labels = np.asarray(subset.var_names.astype(str))[np.asarray(gene_mask, dtype=bool)]
else:
    gene_labels = np.asarray(subset.var_names.astype(str))

for treatment in treatments:
    coef = design.columns.get_loc(treatment)
    res = ep.glm_sc_test(fit, coef=coef)
    results = res["table"] if isinstance(res, dict) else res
    results = results.copy() if hasattr(results, "copy") else pd.DataFrame(results)
    if "genes" not in results.columns:
        if len(gene_labels) != len(results):
            raise ValueError(f"Gene label length ({len(gene_labels)}) does not match results ({len(results)})")
        results.insert(0, "genes", gene_labels)
    out_path = f"{prefix}_{safe_name(treatment)}_results.parquet"
    write_standard_de_parquet(
        results,
        out_path,
        gene_col=None,
        log2fc_col="logFC",
        pvalue_col="PValue",
        padj_col="FDR",
        group=treatment,
        contrast=f"{treatment} vs {reference_condition}",
        stratum=f"celltype={celltype_value}",
    )
    written.append(out_path)

if not written:
    raise ValueError("No edgepython_sc contrasts could be tested")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
