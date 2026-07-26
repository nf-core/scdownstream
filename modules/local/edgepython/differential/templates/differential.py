#!/usr/bin/env python3

import base64
import json
import os
import platform
import re

os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import edgepython as ep
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgePython pseudobulk differential expression")


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


def load_interesting_genes(path_str):
    if not path_str or path_str in ("[]", "null", "None"):
        return set()
    genes = set()
    first_data = True
    with open(path_str) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            field = line.split(",")[0].strip().strip('"').strip("'")
            if not field:
                continue
            if first_data and field.lower() in ("gene", "symbol", "names"):
                first_data = False
                continue
            first_data = False
            genes.add(field.lower())
    return genes


def mqc_parent(method_slug, method_label):
    integration = "${meta.integration}"
    parent_id = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{integration}_{method_slug}")
    return {
        "parent_id": parent_id,
        "parent_name": f"{integration}: {method_label}",
        "parent_description": f"Differential expression volcano plots from {method_label} ({integration} integration).",
    }


def write_volcano(
    df,
    gene_col,
    lfc_col,
    p_col,
    padj_col,
    interesting,
    out_png,
    out_mqc_id,
    section_name,
    description,
    parent,
):
    plot_df = df.copy()
    if gene_col is None:
        if "genes" in plot_df.columns:
            gene_col = "genes"
        else:
            plot_df = plot_df.reset_index()
            gene_col = plot_df.columns[0]
    if gene_col not in plot_df.columns or lfc_col not in plot_df.columns:
        print(f"Warning: missing volcano columns for {out_png}; skipping.")
        return
    use_padj = padj_col is not None and padj_col in plot_df.columns
    p_use = padj_col if use_padj else p_col
    if p_use not in plot_df.columns:
        print(f"Warning: missing p-value column for {out_png}; skipping volcano.")
        return

    plot_df[gene_col] = plot_df[gene_col].astype(str)
    plot_df[lfc_col] = pd.to_numeric(plot_df[lfc_col], errors="coerce")
    plot_df[p_use] = pd.to_numeric(plot_df[p_use], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        plot_df["neglog10"] = -np.log10(plot_df[p_use])
    plot_df["neglog10"] = plot_df["neglog10"].replace([np.inf, -np.inf], np.nan)
    plot_df = plot_df.dropna(subset=[lfc_col, "neglog10"])
    if len(plot_df) < 2:
        print(f"Warning: fewer than 2 plottable points for {out_png}; skipping volcano.")
        return

    if use_padj:
        padj_vals = pd.to_numeric(plot_df[padj_col], errors="coerce")
        significant = (padj_vals < 0.05) & (plot_df[lfc_col].abs() >= 1)
    else:
        significant = (plot_df[p_use] < 0.05) & (plot_df[lfc_col].abs() >= 1)
    plot_df["significant"] = significant.fillna(False)
    plot_df["interesting"] = plot_df[gene_col].str.lower().isin(interesting)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    colours = plot_df["significant"].map({True: "#d62728", False: "#7f7f7f"})
    other = plot_df[~plot_df["interesting"]]
    marked = plot_df[plot_df["interesting"]]
    if not other.empty:
        ax.scatter(
            other[lfc_col], other["neglog10"], c=colours.loc[other.index], marker="o", alpha=0.5, s=16, linewidths=0
        )
    if not marked.empty:
        ax.scatter(
            marked[lfc_col],
            marked["neglog10"],
            c=colours.loc[marked.index],
            marker="^",
            alpha=0.9,
            s=36,
            linewidths=0,
        )
    ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.axvline(-1, color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adjusted p-value)" if use_padj else "-log10(p-value)")
    ax.set_title(section_name)

    label_pool = plot_df[plot_df["significant"]].copy()
    if label_pool.empty:
        label_pool = plot_df.copy()
    label_pool = label_pool.sort_values(p_use, ascending=True)
    prefer = label_pool[label_pool["interesting"]]
    rest = label_pool[~label_pool["interesting"]]
    to_label = pd.concat([prefer, rest]).head(10)
    for _, row in to_label.iterrows():
        ax.annotate(row[gene_col], (row[lfc_col], row["neglog10"]), fontsize=7, alpha=0.9)

    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    with open(out_png, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
    custom_json = {
        "id": out_mqc_id,
        "parent_id": parent["parent_id"],
        "parent_name": parent["parent_name"],
        "parent_description": parent["parent_description"],
        "section_name": section_name,
        "description": description,
        "plot_type": "image",
        "data": image_html,
    }
    with open(f"{out_mqc_id}_mqc.json", "w") as f_json:
        json.dump(custom_json, f_json)


required_cols = ["donor", "condition", "celltype"]
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in pseudobulk adata.obs")

interesting = load_interesting_genes("${interesting_genes}")
volcano_parent = mqc_parent("edgepython", "edgePython")
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
    design = ep.model_matrix(design_formula, metadata)
    y = ep.make_dgelist(counts=counts, samples=metadata)
    y = ep.calc_norm_factors(y)
    y = ep.estimate_disp(y)
    fit = ep.glm_ql_fit(y, design)
    safe_celltype = safe_name(celltype)

    for treatment in treatments:
        coef_name = condition_coef_column(design, treatment)
        if coef_name is None:
            continue
        coef_index = list(design.columns).index(coef_name)
        res = ep.glm_ql_ftest(fit, coef=coef_index)
        top = ep.top_tags(res, n=sub.n_vars)
        out_path = f"{prefix}_{safe_celltype}_{safe_name(treatment)}_results.csv"
        top["table"].to_csv(out_path)
        written.append(out_path)
        stem = out_path.replace("_results.csv", "")
        results_df = pd.read_csv(out_path, index_col=0)
        write_volcano(
            results_df,
            None,
            "logFC",
            "PValue",
            "FDR",
            interesting,
            f"{stem}_volcano.png",
            f"{stem}_volcano",
            f"edgePython volcano: {treatment} vs {reference_condition} (within celltype={celltype})",
            (
                f"edgePython volcano for contrast <code>{treatment}</code> versus <code>{reference_condition}</code> "
                f"within <code>celltype={celltype}</code>."
            ),
            volcano_parent,
        )

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
