#!/usr/bin/env python3

import base64
import json
import platform
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        # An absolute path keeps zarr from parsing colons in Nextflow prefixes as URL schemes
        return ad.read_zarr(Path(path).absolute())
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(Path(f"{stem}.zarr").absolute())


threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

symbol_col = "${symbol_col}"
prefix = "${prefix}"


def read_gene_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


s_genes = read_gene_list("${s_genes}")
g2m_genes = read_gene_list("${g2m_genes}")

adata = read_anndata("${anndata}")

# If gene symbols are in a column rather than the index, temporarily set the
# index so that sc.tl.score_genes_cell_cycle can find them, then restore.
original_index = None
if symbol_col != "index":
    if symbol_col not in adata.var.columns:
        raise ValueError(
            f"symbol_col '{symbol_col}' not found in adata.var. Available columns: {list(adata.var.columns)}"
        )
    original_index = adata.var_names.copy()
    adata.var_names = adata.var[symbol_col].astype(str)
    adata.var_names_make_unique()

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

if original_index is not None:
    adata.var_names = original_index

adata.obs[["S_score", "G2M_score", "phase"]].to_pickle(f"{prefix}.pkl")
write_anndata(adata, prefix, "${anndata_format}")

if adata.obs["phase"].nunique() > 1:
    phase_counts = adata.obs["phase"].value_counts()
    valid_phases = phase_counts[phase_counts > 0].index
    adata_violin = adata[adata.obs["phase"].isin(valid_phases)].copy()
    adata_violin.obs["phase"] = adata_violin.obs["phase"].astype("category").cat.remove_unused_categories()

    sc.pl.violin(adata_violin, ["S_score", "G2M_score"], groupby="phase", show=False)
    fig = plt.gcf()
    plot_path = f"{prefix}_cell_cycle_scores.png"
    plt.savefig(plot_path)

    with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        custom_json = {
            "id": "${prefix}",
            "parent_id": "cell_cycle",
            "parent_name": "Cell cycle scoring",
            "parent_description": "Cell cycle scores assigned before integration.",
            "section_name": "${meta.id}",
            "plot_type": "image",
            "data": image_html,
        }
        json.dump(custom_json, f_json)
else:
    print(f"Warning: only one cell cycle phase detected; skipping violin plot for {prefix}.")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
