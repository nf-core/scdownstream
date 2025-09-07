#!/usr/bin/env python3

import os
import json
import platform
import re
import base64
import pickle
import pathlib

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
    Replaces invalid characters in a filename with underscores.
    """
    # Define a set of characters generally considered invalid in filenames
    # This set covers common invalid characters across Windows and Unix-like systems
    invalid_chars = r'[<>:"/\\|?*\x00-\x1F]' 
    
    # Replace invalid characters with an underscore
    sanitized_filename = re.sub(invalid_chars, '_', filename)
    
    # Remove leading/trailing spaces and dots, which can also cause issues
    sanitized_filename = sanitized_filename.strip(' .')
    
    # Ensure the filename is not empty after sanitization
    if not sanitized_filename:
        sanitized_filename = "untitled" 
        
    return sanitized_filename

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
sample_group_col = "${sample_group_col}"
sample_col = "sample" # always

# read the cluster csv files
cell_groups_csv = "${cluster_csv}"
cell_groups_df = pd.read_csv(cell_groups_csv, index_col=0)
outdir = pathlib.Path(prefix)
os.makedirs(outdir, parents=True, exist_ok=True)

# make sure the annotation file matches the adata
if not set(adata.obs.index) == set(cell_groups_df.index):
    raise ValueError(f"adata has {adata.obs.index.nunique()} cells, cluster_df has {cell_groups_df.index.nunique()} cells")

# ensure they match
cell_cluster_df = cell_groups_df.reindex(adata.obs.index)
cell_group_cols = cell_groups_df.columns

# create new columns for each cluster
for cell_group_col in cell_group_cols:
    if cluster_col in adata.obs.columns:
        raise ValueError(f"Cluster column {cell_group_col} already exists in adata")
    adata.obs[cell_group_col] = cell_groups_df[cell_group_col].astype("str").astype("category")

# We do DE analysis for each cluster column
for cell_group_col in cell_group_cols:
    print(f"Differential analysis for {cell_group_col}")

    # remove groups that have less than 3 cells
    cell_groups_cell_counts = adata.obs[cell_group_col].value_counts()
    cell_groups = cell_groups_cell_counts[cell_groups_cell_counts >= 3].index.astype("str").tolist()
    if len(cell_groups) < 2:
        print(f"    - Skipping {cell_group_col} because it has less than 2 groups with at least 3 cells")
        continue
    
    # set the output directory
    cell_group_outdir = outdir / cell_group_col
    
    #########################################################
    # DE analysis for each group 
    ########################################################
    # for each group, we first compare it to the rest, then we compare it to each other group
    for cell_group in cell_groups:
        print(f"  - Processing cell group {cell_group}:")

        # we make the list to loop over
        other_groups = cell_groups.copy()
        if len(other_groups) > 2:
            # replace the current group with "rest"
            other_groups.remove(cell_group)
            other_groups.append("rest")

        # make the output directory for the one versus the rest DE analysis
        cell_group_outdir = cell_groups_outdir / cell_group

        ##############################################################################
        # Differential analysis #1: Group cells by cell type, then compare between groups
        ##############################################################################
        # for each other group, we compare the current group to it
        for other_group in other_groups:
            print(f"    - Comparing to {other_group}")
            # call the DE analysis
            sc.tl.rank_genes_groups(
                adata,
                groupby = cell_group_col,
                groups = cell_group,
                reference = other_group,
                pts=True)
        
            # get DE results as a dataframe
            rgg_df = sc.get.rank_genes_groups_df(adata, group = None)
        
            # save the rgg_df to a csv file
            os.makedirs(cluster_group_outdir, parents=True, exist_ok=True)
            rgg_df.to_csv(cluster_group_outdir / sanitize_filename(f"{cluster_group}_vs_{other_group}.csv"))

            # plot the rank genes groups
            sc.pl.rank_genes_groups(adata, show=False)
            path = cluster_group_outdir / sanitize_filename(f"{cluster_group}_vs_{other_group}.png")
            plt.savefig(path)
            # close the figure
            plt.close()

        ##############################################################################
        # Differential analysis #2: within each cluster group, compare sample groups
        ##############################################################################
        if sample_group_col == "null":
            continue
        else:
            # make the output directory for the sample group comparisons
            cluster_group_samplegroup_outdir = cluster_group_outdir / "sample_group_comparisons"

            # we take out the current group from the adata
            group_adata = adata[adata.obs[cell_group_col] == cluster_group].copy()
            
            # we take the sample groups that have at least 3 cells
            sample_group_value_counts = group_adata.obs[sample_group_col].astype(str).value_counts()
            sample_group_value_counts = sample_group_value_counts[sample_group_value_counts >= 3]
            sample_groups = sample_group_value_counts.index.astype("str").tolist()
            if len(sample_groups) < 2:
                print(f"    - Skipping sample group comparisons because it has less than 2 groups with at least 3 cells")
                continue

            print(f"    - Comparing between sample groups recorded in {sample_group_col}")
            for sample_group in sample_groups:
                # we make the list to loop over
                other_sample_groups = sample_groups.copy()
                if len(other_sample_groups) > 2:
                    # replace the current sample group with "rest"
                    other_sample_groups.remove(sample_group)
                    other_sample_groups.append("rest")
                
                print(f"      - Processing {sample_group}")
                for other_sample_group in other_sample_groups:
                    print(f"        - Comparing to {other_sample_group}")
                    sc.tl.rank_genes_groups(group_adata, groupby = sample_group_col, groups = sample_group, reference = other_sample_group, pts=True)
                    rgg_df = sc.get.rank_genes_groups_df(group_adata, group = None)
                    rgg_df.to_csv(cluster_group_samplegroup_outdir / sanitize_filename(f"{sample_group}_vs_{other_sample_group}.csv"))

                    sc.pl.rank_genes_groups(group_adata, show=False)
                    path = cluster_group_samplegroup_outdir / sanitize_filename(f"{sample_group}_vs_{other_sample_group}.png")
                    plt.savefig(path)

    #########################################################
    # DE analysis for each sample group 
    ########################################################
    if sample_group_col == "null":
        continue
    else:
        # we take the sample groups that have at least 3 cells
        sample_group_value_counts = adata.obs[sample_group_col].astype(str).value_counts()
        sample_group_value_counts = sample_group_value_counts[sample_group_value_counts >= 3]
        sample_groups = sample_group_value_counts.index.astype("str").tolist()
        if len(sample_groups) < 2:
            print(f"    - Skipping sample group comparisons because it has less than 2 groups with at least 3 cells")
            continue

        # for each sample group, we first compare it to the rest, then we compare it to each other sample group
        for sample_group in sample_groups:
            print(f"  - Processing sample group {sample_group}:")

            # we make the list to loop over
            other_groups = sample_groups.copy()
            if len(other_groups) > 2:
                # replace the current group with "rest"
                other_groups.remove(sample_group)
                other_groups.append("rest")

            # make the output directory for the one versus the rest DE analysis
            sample_group_outdir = outdir / sample_group

            ##############################################################################
            # Differential analysis #1: Group cells by sample group, then compare between groups
            ##############################################################################
            # for each other group, we compare the current group to it
            for other_group in other_groups:
                print(f"    - Comparing to {other_group}")
                # call the DE analysis
                sc.tl.rank_genes_groups(
                    adata,
                    groupby = cell_group_col,
                    groups = sample_group,
                    reference = other_group,
                    pts=True)
            
                # get DE results as a dataframe
                rgg_df = sc.get.rank_genes_groups_df(adata, group = None)
            
                # save the rgg_df to a csv file
                os.makedirs(sample_group_outdir, parents=True, exist_ok=True)
                rgg_df.to_csv(sample_group_outdir / sanitize_filename(f"{sample_group}_vs_{other_group}.csv"))

                # plot the rank genes groups
                sc.pl.rank_genes_groups(adata, show=False)
                path = sample_group_outdir / sanitize_filename(f"{sample_group}_vs_{other_group}.png")
                plt.savefig(path)
                # close the figure
                plt.close()

            ##############################################################################
            # Differential analysis #2: within each sample group, compare cluster groups
            ##############################################################################

            # make the output directory for the sample group comparisons
            sample_group_clustergroup_outdir = sample_group_outdir / "cluster_group_comparisons"

            # we take out the current group from the adata
            group_adata = adata[adata.obs[cell_group_col] == sample_group].copy()
            
            # we take the cluster groups that have at least 3 cells
            sample_group_value_counts = group_adata.obs[cell_group_col].astype(str).value_counts()
            sample_group_value_counts = sample_group_value_counts[sample_group_value_counts >= 3]
            cluster_groups = sample_group_value_counts.index.astype("str").tolist()
            if len(cluster_groups) < 2:
                print(f"    - Skipping sample group comparisons because it has less than 2 groups with at least 3 cells")
                continue

            print(f"    - Comparing between sample groups recorded in {sample_group_col}")
            for sample_group in sample_groups:
                # we make the list to loop over
                other_sample_groups = sample_groups.copy()
                if len(other_sample_groups) > 2:
                    # replace the current sample group with "rest"
                    other_sample_groups.remove(sample_group)
                    other_sample_groups.append("rest")
                
                print(f"      - Processing {sample_group}")
                for other_sample_group in other_sample_groups:
                    print(f"        - Comparing to {other_sample_group}")
                    sc.tl.rank_genes_groups(group_adata, groupby = sample_group_col, groups = sample_group, reference = other_sample_group, pts=True)
                    rgg_df = sc.get.rank_genes_groups_df(group_adata, group = None)
                    rgg_df.to_csv(sample_group_samplegroup_outdir / sanitize_filename(f"{sample_group}_vs_{other_sample_group}.csv"))

                    sc.pl.rank_genes_groups(group_adata, show=False)
                    path = sample_group_samplegroup_outdir / sanitize_filename(f"{sample_group}_vs_{other_sample_group}.png")
                    plt.savefig(path)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)


