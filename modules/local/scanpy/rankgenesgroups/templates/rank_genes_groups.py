#!/usr/bin/env python3

import os
import json
import platform
import base64
import pickle

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

def sanitize_filename(filename):
    """
    Sanitize filename by replacing special characters with underscores.
    Keeps alphanumeric characters and underscores, replaces everything else.
    """
    if not filename:
        return "unknown"
    
    # Replace spaces and hyphens with underscores first
    sanitized = filename.replace(" ", "_").replace("-", "_")
    
    # Keep only alphanumeric characters and underscores
    sanitized = "".join([c if c.isalnum() or c == "_" else "_" for c in sanitized])
    
    # Remove multiple consecutive underscores
    while "__" in sanitized:
        sanitized = sanitized.replace("__", "_")
    
    # Remove leading/trailing underscores
    sanitized = sanitized.strip("_")
    
    # If empty after sanitization, return default
    if not sanitized:
        return "unknown"
    
    return sanitized

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
sample_group_col = "${sample_group_col}"
sample_col = "sample" # always

# read the cluster csv files
cluster_csv = "${cluster_csv}"
cluster_df = pd.read_csv(cluster_csv, index_col=0)
cluster_col = cluster_df.columns[0]
outdir = prefix
os.makedirs(outdir, exist_ok=True)

# check if adata and cluster_df have the same index
if not adata.obs.index.sort_values().equals(cluster_df.index.sort_values()):
    print(f"warning: adata has {adata.obs.index.nunique()} cells, cluster_df has {cluster_df.index.nunique()} cells")

# ensure they match
cluster_df = cluster_df[cluster_df.index.isin(adata.obs.index)]
adata = adata[adata.obs.index.isin(cluster_df.index)]
cluster_df = cluster_df.reindex(adata.obs.index)
adata.obs[cluster_col] = cluster_df[cluster_col].astype("str").astype("category")

# Differential analysis #1: pair-wise cell type comparisons
print("Differential analysis #1: pair-wise cell type comparisons")

# keep groups that have at least 3 cells
cluster_cell_counts = adata.obs[cluster_col].value_counts()
cluster_groups = cluster_cell_counts[cluster_cell_counts >= 3].index.astype("str").tolist()

sc.tl.rank_genes_groups(adata, groupby = cluster_col, groups = cluster_groups, pts=True)
rgg_dict = adata.uns["rank_genes_groups"]
rgg_df = sc.get.rank_genes_groups_df(adata, group = None)

# save the rgg_df to a csv file
rgg_df.to_csv(f"{outdir}/pairwise_comparisons.csv")

# plot the rank genes groups
sc.pl.rank_genes_groups(adata, show=False)
path = f"{outdir}/pairwise_comparisons.png"
plt.savefig(path)

# Differential analysis #2: within each cluster group, compare samples
print("Differential analysis #2: within each cluster group, compare samples")
cluster_sample_col = "cluster_sample"
os.makedirs(f"{outdir}/{cluster_sample_col}", exist_ok=True)
adata.obs[cluster_sample_col] = (adata.obs[cluster_col].astype(str) + "_" + adata.obs[sample_col].astype(str)).astype("category")

cluster_sample_value_counts_all = adata.obs[cluster_sample_col].astype(str).value_counts()

for cluster in cluster_groups:
    cluster_sanitized = sanitize_filename(str(cluster))
    os.makedirs(f"{outdir}/{cluster_sanitized}", exist_ok=True)
    # we define the `groups` argument to include only the current cluster
    cluster_sample_groups = cluster_sample_value_counts_all[cluster_sample_value_counts_all.index.str.startswith(str(cluster)+'_')]
    cluster_sample_groups = cluster_sample_groups[cluster_sample_groups >= 3].index.astype("str").tolist()
    
    if len(cluster_sample_groups) >= 2:
        sc.tl.rank_genes_groups(adata, groupby = cluster_sample_col, groups = cluster_sample_groups, pts=True)
        rgg_df = sc.get.rank_genes_groups_df(adata, group = None)

        rgg_df.to_csv(f"{outdir}/{cluster_sample_col}/{cluster_sanitized}.csv")

        sc.pl.rank_genes_groups(adata, show=False)
        path = f"{outdir}/{cluster_sample_col}/{cluster_sanitized}.png"
        plt.savefig(path)

# Differential analysis #3: within each cluster group, compare sample groups
if sample_group_col != "null":
    print("Differential analysis #3: within each cluster group, compare sample groups")
    cluster_sample_group_col = "cluster_samplegroup"
    os.makedirs(f"{outdir}/{cluster_sample_group_col}", exist_ok=True)
    adata.obs[cluster_sample_group_col] = (adata.obs[cluster_col].astype(str) + "_" + adata.obs[sample_group_col].astype(str)).astype("category")
    
    cluster_sample_group_groups_value_counts_all = adata.obs[cluster_sample_group_col].astype(str).value_counts()

    for cluster_sample_group in cluster_sample_group_groups_value_counts_all.index:
        cluster_sample_group_sanitized = sanitize_filename(str(cluster_sample_group))
        os.makedirs(f"{outdir}/{cluster_sample_group_sanitized}", exist_ok=True)
        cluster_sample_group_groups = cluster_sample_group_groups_value_counts_all[cluster_sample_group_groups_value_counts_all.index.str.startswith(str(cluster_sample_group)+'_')]
        cluster_sample_group_groups = cluster_sample_group_groups[cluster_sample_group_groups >= 3].index.astype("str").tolist()
        if len(cluster_sample_group_groups) >= 2:
            sc.tl.rank_genes_groups(adata, groupby = cluster_sample_group_col, groups = cluster_sample_group_groups, pts=True)
            rgg_df = sc.get.rank_genes_groups_df(adata, group = None)
            rgg_df.to_csv(f"{outdir}/{cluster_sample_group_col}/{cluster_sample_group_sanitized}.csv")

            sc.pl.rank_genes_groups(adata, show=False)
            path = f"{outdir}/{cluster_sample_group_col}/{cluster_sample_group_sanitized}.png"
            plt.savefig(path)
else:
    print("Skipping sample group comparison as sample_group_col is 'none'")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
