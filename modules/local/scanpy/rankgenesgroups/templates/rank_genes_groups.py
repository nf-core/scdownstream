#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import base64
import json
import pickle
import platform
import re

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
method = "${method}"
rank_key = "${rank_key}"

filter_col = "${filter_col ?: ''}"
filter_val = "${filter_val ?: ''}"

meta_id = "${meta.id}"
obs_key = "${obs_key}"

adata.obs[obs_key] = adata.obs[obs_key].astype(str)

if filter_col and filter_val:
    adata = adata[adata.obs[filter_col] == filter_val].copy()

kwargs = {"groupby": obs_key, "method": method, "pts": True, "key_added": rank_key}
filtered_rank_key = f"{rank_key}_filtered"


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


def prepare_volcano_df(df, gene_col, lfc_col, p_col, padj_col, interesting):
    plot_df = df.copy()
    if gene_col is None:
        plot_df = plot_df.reset_index()
        gene_col = plot_df.columns[0]
    if gene_col not in plot_df.columns or lfc_col not in plot_df.columns:
        return None, None, None, None
    use_padj = padj_col is not None and padj_col in plot_df.columns
    p_use = padj_col if use_padj else p_col
    if p_use not in plot_df.columns:
        return None, None, None, None

    plot_df[gene_col] = plot_df[gene_col].astype(str)
    plot_df[lfc_col] = pd.to_numeric(plot_df[lfc_col], errors="coerce")
    plot_df[p_use] = pd.to_numeric(plot_df[p_use], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        plot_df["neglog10"] = -np.log10(plot_df[p_use])
    plot_df["neglog10"] = plot_df["neglog10"].replace([np.inf, -np.inf], np.nan)
    plot_df = plot_df.dropna(subset=[lfc_col, "neglog10"])
    if len(plot_df) < 2:
        return None, None, None, None

    if use_padj:
        padj_vals = pd.to_numeric(plot_df[padj_col], errors="coerce")
        significant = (padj_vals < 0.05) & (plot_df[lfc_col].abs() >= 1)
    else:
        significant = (plot_df[p_use] < 0.05) & (plot_df[lfc_col].abs() >= 1)
    plot_df["significant"] = significant.fillna(False)
    plot_df["interesting"] = plot_df[gene_col].str.lower().isin(interesting)
    return plot_df, gene_col, p_use, use_padj


def draw_volcano_ax(ax, plot_df, gene_col, lfc_col, p_use, use_padj, title, max_labels=3):
    colours = plot_df["significant"].map({True: "#d62728", False: "#7f7f7f"})
    other = plot_df[~plot_df["interesting"]]
    marked = plot_df[plot_df["interesting"]]
    if not other.empty:
        ax.scatter(
            other[lfc_col], other["neglog10"], c=colours.loc[other.index], marker="o", alpha=0.5, s=10, linewidths=0
        )
    if not marked.empty:
        ax.scatter(
            marked[lfc_col],
            marked["neglog10"],
            c=colours.loc[marked.index],
            marker="^",
            alpha=0.9,
            s=28,
            linewidths=0,
        )
    ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.axvline(-1, color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.axvline(1, color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.set_xlabel("log2FC", fontsize=8)
    ax.set_ylabel("-log10(padj)" if use_padj else "-log10(p)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)

    label_pool = plot_df[plot_df["significant"]].copy()
    if label_pool.empty:
        label_pool = plot_df.copy()
    label_pool = label_pool.sort_values(p_use, ascending=True)
    prefer = label_pool[label_pool["interesting"]]
    rest = label_pool[~label_pool["interesting"]]
    to_label = pd.concat([prefer, rest]).head(max_labels)
    for _, row in to_label.iterrows():
        ax.annotate(row[gene_col], (row[lfc_col], row["neglog10"]), fontsize=6, alpha=0.9)


def write_volcano_grid(
    group_frames,
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
    panels = []
    for group_name, gdf in group_frames:
        plot_df, resolved_gene_col, p_use, use_padj = prepare_volcano_df(
            gdf, gene_col, lfc_col, p_col, padj_col, interesting
        )
        if plot_df is None:
            continue
        panels.append((group_name, plot_df, resolved_gene_col, p_use, use_padj))

    if not panels:
        print(f"Warning: no plottable volcano panels for {out_png}; skipping.")
        return

    n = len(panels)
    ncols = min(4, max(1, int(np.ceil(np.sqrt(n)))))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.6 * ncols, 3.0 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    for idx, (group_name, plot_df, resolved_gene_col, p_use, use_padj) in enumerate(panels):
        row, col = divmod(idx, ncols)
        draw_volcano_ax(
            axes[row][col],
            plot_df,
            resolved_gene_col,
            lfc_col,
            p_use,
            use_padj,
            f"{group_name} vs rest",
        )
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    fig.suptitle(section_name, fontsize=11)
    plt.savefig(out_png, bbox_inches="tight", dpi=150)
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


# Check value counts for each group
value_counts = adata.obs[obs_key].value_counts()
# Filter out groups with less than 2 samples (scanpy requirement)
valid_groups = value_counts[value_counts >= 2].index.tolist()
invalid_groups = value_counts[value_counts < 2].index.tolist()

if len(invalid_groups) > 0:
    print(f"Warning: Excluding groups with < 2 samples: {', '.join(map(str, invalid_groups))}")

# Only proceed if we have at least 2 groups with >= 2 samples each
if len(valid_groups) >= 2:
    # Filter adata to only include valid groups
    adata = adata[adata.obs[obs_key].isin(valid_groups)].copy()

    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, **kwargs)
    sc.tl.filter_rank_genes_groups(
        adata,
        key=rank_key,
        key_added=filtered_rank_key,
        min_in_group_fraction=0.2,
        max_out_group_fraction=0.2,
    )

    interesting = load_interesting_genes("${interesting_genes}")
    method_label = f"Scanpy {method}"
    method_slug = re.sub(r"[^A-Za-z0-9._-]+", "_", f"scanpy_{method}")
    volcano_parent = mqc_parent(method_slug, method_label)
    try:
        full_df = sc.get.rank_genes_groups_df(adata, group=None, key=rank_key)
        volcano_section = f"{method_label} volcanoes: {obs_key} vs rest"
        volcano_description = (
            f"{method_label} volcano panels for each <code>{obs_key}</code> group versus the rest of the cells."
        )
        if filter_col and filter_val:
            volcano_section = f"{method_label} volcanoes: {obs_key} vs rest (within {filter_col}={filter_val})"
            volcano_description = (
                f"{method_label} volcano panels for each <code>{obs_key}</code> group versus the rest of the cells "
                f"within <code>{filter_col}={filter_val}</code>."
            )
        write_volcano_grid(
            list(full_df.groupby("group", sort=True)),
            "names",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            interesting,
            f"{prefix}_volcano.png",
            f"{prefix}_volcano",
            volcano_section,
            volcano_description,
            volcano_parent,
        )
    except Exception as exc:
        print(f"Warning: skipping volcano plots: {exc}")

    marker_df = sc.get.rank_genes_groups_df(adata, group=None, key=filtered_rank_key)
    marker_df = marker_df[marker_df["names"].notna()].copy()

    if marker_df.empty:
        print(f"Warning: no genes passed filter for {obs_key}; skipping plots and H5AD output.")
        adata.uns.pop(filtered_rank_key, None)
        adata.uns.pop(rank_key, None)
    else:
        marker_df.to_csv(f"{prefix}_markers.csv", index=False)

        # Store filtered markers under rank_key for downstream consumers such as CyteType.
        rgg_dict = dict(adata.uns.pop(filtered_rank_key))
        adata.uns.pop(rank_key, None)
        # Scanpy marks filtered-out genes as missing names; replace them so H5AD serialisation works.
        names_df = pd.DataFrame(rgg_dict["names"]).fillna("").astype(str)
        rgg_dict["names"] = names_df.to_records(index=False)
        adata.uns[rank_key] = rgg_dict

        pickle.dump(rgg_dict, open(f"{prefix}.pkl", "wb"))
        adata.write_h5ad(f"{prefix}.h5ad")

        # Plot
        sc.pl.rank_genes_groups(adata, key=rank_key, show=False)
        path = f"{prefix}.png"
        plt.savefig(path, bbox_inches="tight")

        sc.pl.rank_genes_groups_dotplot(adata, key=rank_key, show=False)
        dotplot_path = f"{prefix}_dotplot.png"
        plt.savefig(dotplot_path, bbox_inches="tight")

        # Build section name with filter and obs_key information
        if filter_col and filter_val:
            section_name = f"Characteristic genes (grouped by: {obs_key}, filtered: {filter_col}={filter_val})"
            description = f"Characteristic genes, grouped by <code>{obs_key}</code>, filtered to <code>{filter_col}={filter_val}</code>."
        else:
            section_name = f"Characteristic genes (grouped by: {obs_key})"
            description = f"Characteristic genes, grouped by <code>{obs_key}</code>."

        def write_mqc_plot(plot_path, plot_id, plot_label):
            with open(plot_path, "rb") as f_plot:
                image_string = base64.b64encode(f_plot.read()).decode("utf-8")
            image_html = (
                f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
            )
            custom_json = {
                "id": plot_id,
                "parent_id": "${meta.integration}",
                "parent_name": "${meta.integration}",
                "parent_description": "Results of the ${meta.integration} integration.",
                "section_name": f"{section_name} ({plot_label})",
                "description": f"{description} {plot_label.capitalize()}.",
                "plot_type": "image",
                "data": image_html,
            }
            with open(f"{plot_id}_mqc.json", "w") as f_json:
                json.dump(custom_json, f_json)

        write_mqc_plot(path, "${prefix}", "rank plot")
        write_mqc_plot(dotplot_path, "${prefix}_dotplot", "dot plot")
else:
    if len(valid_groups) == 0:
        print("Skipping rank_genes_groups computation: no groups have >= 2 samples.")
    elif len(valid_groups) == 1:
        print(f"Skipping rank_genes_groups computation: only one group has >= 2 samples (group: {valid_groups[0]}).")
    else:
        print("Skipping rank_genes_groups computation: less than 2 valid groups remaining after filtering.")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
