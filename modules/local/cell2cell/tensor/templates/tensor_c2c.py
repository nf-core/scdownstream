#!/usr/bin/env python3

import base64
import json
import os
import platform
import re

import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"
os.environ["OMP_NUM_THREADS"] = "${task.cpus}"

import matplotlib

matplotlib.use("Agg")
import cell2cell as c2c
import decoupler as dc
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

prefix = "${prefix}"
integration = "${integration}"
context_key = "${context_key}"
species = "${species}"
rank_raw = "${rank}"
seed = int("${seed}")


def _optional_int(value):
    return None if value in ("", "null", "None") else int(value)


def _slug(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def _write_mqc(png_path, section_name, plot_id=None):
    if plot_id is None:
        plot_id = os.path.basename(png_path)
        if plot_id.endswith(".png"):
            plot_id = plot_id[:-4]
    with open(png_path, "rb") as f_plot, open(f"{plot_id}_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        json.dump(
            {
                "id": plot_id,
                "parent_id": integration,
                "parent_name": integration,
                "parent_description": f"Results of the {integration} integration.",
                "section_name": section_name,
                "plot_type": "image",
                "data": image_html,
            },
            f_json,
        )


rank = _optional_int(rank_raw)

liana_res = pd.read_csv("${liana_bysample}")
contexts = pd.read_csv("${contexts}", sep="\t")

if context_key not in liana_res.columns:
    raise SystemExit(
        f"Context column '{context_key}' missing from by-sample results; columns: {list(liana_res.columns)}"
    )
if context_key not in contexts.columns:
    raise SystemExit(f"Context column '{context_key}' missing from contexts table; columns: {list(contexts.columns)}")
if liana_res[context_key].nunique() < 2:
    raise SystemExit(
        f"At least two contexts are required for tensor-cell2cell; found {liana_res[context_key].nunique()}."
    )
if "magnitude_rank" not in liana_res.columns:
    raise SystemExit(f"Column 'magnitude_rank' missing from by-sample results; columns: {list(liana_res.columns)}")

context_dict = None
if "condition" in contexts.columns:
    grouped = contexts.groupby(context_key, observed=True)["condition"].nunique()
    multi = grouped[grouped > 1]
    if not multi.empty:
        print(
            "Warning: omitting condition metadata because some contexts "
            f"map to multiple conditions: {multi.index.tolist()}"
        )
    else:
        context_dict = {str(row[context_key]): str(row["condition"]) for _, row in contexts.iterrows()}

tensor = li.multi.to_tensor_c2c(
    liana_res=liana_res,
    sample_key=context_key,
    score_key="magnitude_rank",
    how="outer_cells",
)

metadata_dicts = [context_dict, None, None, None]
tensor_meta = c2c.tensor.generate_tensor_metadata(
    interaction_tensor=tensor,
    metadata_dicts=metadata_dicts,
    fill_with_order_elements=True,
)

tensor = c2c.analysis.run_tensor_cell2cell_pipeline(
    tensor,
    tensor_meta,
    copy_tensor=True,
    rank=rank,
    tf_optimization="regular",
    random_state=seed,
    device="cpu",
    elbow_metric="error",
    smooth_elbow=False,
    upper_rank=25,
    tf_init="random",
    tf_svd="numpy_svd",
    output_fig=False,
)

factors = tensor.factors
factor_names = list(factors["Sender Cells"].columns)

dimension_files = {
    "Contexts": f"{prefix}_loadings_contexts.csv",
    "Ligand-Receptor Pairs": f"{prefix}_loadings_ligand_receptor_pairs.csv",
    "Sender Cells": f"{prefix}_loadings_sender_cells.csv",
    "Receiver Cells": f"{prefix}_loadings_receiver_cells.csv",
}

for dim_name, filename in dimension_files.items():
    if dim_name not in factors:
        raise SystemExit(f"Missing factor dimension '{dim_name}' in tensor factors")
    factors[dim_name].to_csv(filename)

# Factor overview (R plot_c2c_overview)
overview_png = f"{prefix}_tensor_factors.png"
c2c.plotting.tensor_factors_plot(
    interaction_tensor=tensor,
    metadata=tensor_meta,
    sample_col="Element",
    group_col="Category",
    meta_cmaps=["viridis", "Dark2_r", "tab20", "tab20"],
    fontsize=10,
    filename=overview_png,
)
plt.close("all")
_write_mqc(overview_png, "${meta.id} Cell2cell factor overview")

# LR loadings heatmap (R plot_lr_heatmap)
lr_png = f"{prefix}_loadings_lr_clustermap.png"
lr_loadings = factors["Ligand-Receptor Pairs"]
lr_threshold = 0.1
if not (lr_loadings > lr_threshold).any(axis=1).any():
    # Keep the strongest LRs when all loadings fall below the default threshold
    top_n = min(25, lr_loadings.shape[0])
    keep = lr_loadings.max(axis=1).nlargest(top_n).index
    lr_loadings = lr_loadings.loc[keep]
    lr_threshold = 0.0
c2c.plotting.loading_clustermap(
    loadings=lr_loadings,
    loading_threshold=lr_threshold,
    use_zscore=False,
    figsize=(28, 8),
    cmap="RdBu_r",
    cbar_label="LR Loadings",
    row_cluster=False,
    filename=lr_png,
)
plt.close("all")
_write_mqc(lr_png, "${meta.id} Cell2cell LR loadings")

# Context loadings heatmap (R plot_context_heat)
contexts_png = f"{prefix}_loadings_contexts_clustermap.png"
c2c.plotting.loading_clustermap(
    loadings=factors["Contexts"],
    use_zscore=False,
    figsize=(16, 6),
    cmap="RdBu_r",
    cbar_label="Context Loadings",
    filename=contexts_png,
)
plt.close("all")
_write_mqc(contexts_png, "${meta.id} Cell2cell context loadings")

# Context boxplots (R plot_context_boxplot); needs unique condition labels
if context_dict and len(set(context_dict.values())) >= 2:
    box_png = f"{prefix}_context_boxplots.png"
    try:
        c2c.plotting.context_boxplot(
            context_loadings=factors["Contexts"],
            metadict=context_dict,
            nrows=2,
            figsize=(12, 6),
            filename=box_png,
        )
        plt.close("all")
        _write_mqc(box_png, "${meta.id} Cell2cell context boxplots")
    except Exception as exc:
        print(f"Warning: skipped context boxplots: {exc}")
        plt.close("all")

# PROGENy pathway enrichment dotplot (R vignette enrichment dotplot)
pathway_png = f"{prefix}_pathway_enrichment_dotplot.png"
try:
    lr_pairs = factors["Ligand-Receptor Pairs"]
    example_lr = str(lr_pairs.index[0])
    lr_sep = "^" if "^" in example_lr else ("_" if "_" in example_lr else "^")
    organism_key = species.strip().lower()
    if organism_key in ("homo_sapiens", "hsapiens", "hs", ""):
        organism_key = "human"
    elif organism_key in ("mus_musculus", "mmusculus", "mm"):
        organism_key = "mouse"

    net = dc.op.progeny(organism=organism_key, top=5000)
    resource = li.resource.select_resource("consensus")
    lr_progeny = li.resource.generate_lr_geneset(resource, net, lr_sep=lr_sep)
    if "interaction" in lr_progeny.columns and "target" not in lr_progeny.columns:
        lr_progeny = lr_progeny.rename(columns={"interaction": "target"})

    estimate, pvals = dc.mt.ulm(lr_pairs.T, lr_progeny, tmin=5, raw=False)
    enrich_df = (
        estimate.melt(ignore_index=False, var_name="pathway", value_name="score")
        .rename_axis("factor")
        .reset_index()
        .merge(
            pvals.melt(ignore_index=False, var_name="pathway", value_name="pval").rename_axis("factor").reset_index(),
            on=["factor", "pathway"],
        )
    )
    enrich_df["size"] = -np.log10(enrich_df["pval"].clip(lower=1e-36))
    enrich_df["significant"] = enrich_df["pval"] <= 0.05
    enrich_df.to_csv(f"{prefix}_pathway_enrichment.csv", index=False)

    pathways = sorted(enrich_df["pathway"].unique())
    factors_ord = list(estimate.index)
    enrich_df["x"] = enrich_df["pathway"].map({p: i for i, p in enumerate(pathways)})
    enrich_df["y"] = enrich_df["factor"].map({f: i for i, f in enumerate(factors_ord)})

    sizes = enrich_df["size"]
    size_min, size_max = float(sizes.min()), float(sizes.max())
    if size_max <= size_min:
        point_sizes = np.full(len(enrich_df), 40.0)
    else:
        point_sizes = 20.0 + 180.0 * (sizes - size_min) / (size_max - size_min)

    fig_w = max(8.0, len(pathways) * 0.7)
    fig_h = max(3.5, len(factors_ord) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    scatter = ax.scatter(
        enrich_df["x"],
        enrich_df["y"],
        c=enrich_df["score"],
        s=point_sizes,
        cmap="Reds",
        edgecolors="black",
        linewidths=np.where(enrich_df["significant"], 0.8, 0.3),
        alpha=0.9,
    )
    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(pathways, rotation=90, ha="center")
    ax.set_yticks(range(len(factors_ord)))
    ax.set_yticklabels(factors_ord)
    ax.set_xlabel("Pathway")
    ax.set_ylabel("Factor")
    ax.set_title("PROGENy pathway enrichment")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Activity")
    fig.tight_layout()
    fig.savefig(pathway_png, dpi=300, bbox_inches="tight")
    plt.close("all")
    _write_mqc(pathway_png, "${meta.id} Cell2cell pathway enrichment")
except Exception as exc:
    print(f"Warning: skipped pathway enrichment dotplot: {exc}")
    plt.close("all")

# Sender x receiver loadings-product heatmaps (R plot_c2c_cells)
for factor_name in factor_names:
    loading_product = c2c.analysis.tensor_downstream.get_joint_loadings(
        factors,
        dim1="Sender Cells",
        dim2="Receiver Cells",
        factor=factor_name,
    )
    factor_slug = _slug(factor_name)
    outfile = f"{prefix}_{factor_slug}_loadings_product.png"
    c2c.plotting.loading_clustermap(
        loading_product.T,
        use_zscore=False,
        figsize=(8, 8),
        cmap="RdBu_r",
        cbar_label="Loading Product",
        filename=outfile,
    )
    plt.close("all")
    _write_mqc(
        outfile,
        "${meta.id} Cell2cell {factor_name}",
        plot_id=f"{prefix}_{factor_slug}",
    )

c2c.io.export_variable_with_pickle(tensor, f"{prefix}_tensor.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "liana": li.__version__,
        "cell2cell": c2c.__version__,
        "decoupler": dc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
