#!/usr/bin/env python3

# In this script, we do the following:
# For each cell group column, 
# 1. we do the following fr each cell group recorded in the cell group column:
#    - per-group DE vs each other group and vs rest
#    - if sample groups are provided, we take out the cells with the cell group, and do per-sample group DE vs each other sample group and vs rest
# 2. If sample groups are provided, we do the following for each sample group:
#    - subset the cells with the sample group
#    - then compare cell groups within that subset (for this cell_group_col)

import os
import platform
import re
from pathlib import Path

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


def sanitize_filename(filename: str) -> str:
    """
    Replaces invalid characters in a filename with underscores.
    """
    # Build the pattern programmatically to avoid relying on backslash escapes
    # that could be altered by Nextflow template processing.
    invalid_specials = re.escape('<>:"/' + chr(92) + '|?*')
    control_chars = ''.join(chr(c) for c in range(0x00, 0x20))
    invalid_chars_pattern = f'[{invalid_specials}{control_chars}]'
    sanitized_filename = re.sub(invalid_chars_pattern, "_", str(filename))
    sanitized_filename = sanitized_filename.strip(" .")
    if not sanitized_filename:
        sanitized_filename = "untitled"
    return sanitized_filename


def ensure_categorical_str(adata_obj: sc.AnnData, column: str) -> None:
    """
    Ensures the specified column in the AnnData object is a categorical string.
    """
    adata_obj.obs[column] = adata_obj.obs[column].astype(str).astype("category")


def valid_groups(adata_obj: sc.AnnData, column: str, min_cells: int = 3):
    """
    Returns the valid groups in the specified column of the AnnData object that have at least min_cells cells
    """
    vc = adata_obj.obs[column].astype(str).value_counts()
    return vc[vc >= min_cells].index.astype("str").tolist()


def run_and_save_de(adata_obj: sc.AnnData, groupby: str, group: str, reference, out_dir: Path) -> None:
    """
    Runs differential expression analysis for the specified group and reference group, and saves the results to the specified output directory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    # Run differential expression analysis
    sc.tl.rank_genes_groups(
        adata_obj,
        use_raw=False,
        groupby=groupby,
        groups=[group],
        reference=reference,
        pts=True,
    )
    # Get the results of the differential expression analysis
    rgg_df = sc.get.rank_genes_groups_df(adata_obj, group=None)
    # Create a standardized filename for the output files
    ref_name = reference if isinstance(reference, str) else str(reference)
    basename = f"{sanitize_filename(group)}_vs_{sanitize_filename(ref_name)}"
    rgg_df.to_csv(out_dir / f"{basename}.csv", index=False)
    # Save the results of the differential expression analysis to a png file
    sc.pl.rank_genes_groups(adata_obj, show=False)
    plt.savefig(out_dir / f"{basename}.png")
    plt.close()

adata = sc.read_h5ad("${h5ad}")
sample_group_col = "${sample_group_col}"

cell_groups_csv = "${cluster_csv}"
cell_groups_df = pd.read_csv(cell_groups_csv, index_col=0)
outdir = Path("${prefix}")
outdir.mkdir(parents=True, exist_ok=True)

# Align order
cell_groups_df = cell_groups_df[cell_groups_df.index.isin(adata.obs.index)]
adata = adata[adata.obs.index.isin(cell_groups_df.index)]
cell_groups_df = cell_groups_df.reindex(adata.obs.index)
cell_group_cols = list(cell_groups_df.columns)

# ensure there are cells left
if adata.n_obs == 0 or cell_groups_df.shape[0] == 0:
    ValueError(f"No cells left after aligning adata and cluster_csv")

# Add grouping columns to adata.obs
for cell_group_col in cell_group_cols:
    if cell_group_col not in adata.obs.columns:
        adata.obs[cell_group_col] = cell_groups_df[cell_group_col].astype("str").astype("category")

has_sample_groups = sample_group_col not in (None, "", "null")
if has_sample_groups:
    if sample_group_col not in adata.obs.columns:
        raise ValueError(f"sample_group_col '{sample_group_col}' not found in adata.obs")
    # if all cells have the same sample group, we skip the sample group comparisons
    if adata.obs[sample_group_col].nunique() == 1:
        has_sample_groups = False
        print(f"All cells have the same sample group")
else:
    print(f"No sample group column provided")
if not has_sample_groups:
    print(f"Skipping sample group comparisons.")

for cell_group_col in cell_group_cols:
    print(f"Differential analysis for cell group column: {cell_group_col}")
    ensure_categorical_str(adata, cell_group_col)
    col_outdir = outdir / sanitize_filename(cell_group_col)

    # keep only the groups that have at least 3 cells
    groups = valid_groups(adata, cell_group_col, min_cells=3)
    # if there are fewer than 2 groups with at least 3 cells, we skip the differential analysis
    if len(groups) < 2:
        print(f"\t- Skipping {cell_group_col}: fewer than 2 groups with ≥3 cells")
        continue

    for group in groups:
        # ------------------------------------------------------------
        # 1) For each cell group column: per-group DE vs each other group and vs rest
        # ------------------------------------------------------------
        print(f"\t- Processing cell group '{group}'")
        group_outdir = col_outdir / "cell_groups" / sanitize_filename(group)
        # Pairwise comparisons
        for other in [g for g in groups if g != group]:
            print(f"\t\t- {group} vs {other}")
            run_and_save_de(adata, cell_group_col, group, other, group_outdir)
        # Versus rest (only meaningful if >2 groups)
        if len(groups) > 2:
            print(f"\t\t- {group} vs rest")
            run_and_save_de(adata, cell_group_col, group, "rest", group_outdir)

        # ------------------------------------------------------------
        # 2) For each cell group column: per-sample group DE vs each other sample group and vs rest
        # ------------------------------------------------------------
        if has_sample_groups:
            # 2) Within each cell group, compare sample groups (if provided)
            subset = adata[adata.obs[cell_group_col].astype(str) == group].copy()
            ensure_categorical_str(subset, cell_group_col)
            ensure_categorical_str(subset, sample_group_col)

            sample_groups = valid_groups(subset, sample_group_col, min_cells=3)
            if len(sample_groups) < 2:
                print("\t\t- Skipping within-group sample comparisons: fewer than 2 sample groups with ≥3 cells")
                continue
            for sg in sample_groups:
                print(f"\t\t- Sample group '{sg}' within cell group '{group}'")
                sg_outdir = group_outdir / sanitize_filename(sg)
                for other in [x for x in sample_groups if x != sg]:
                    print(f"\t\t\t- {sg} vs {other}")
                    run_and_save_de(subset, sample_group_col, sg, other, sg_outdir)
                if len(sample_groups) > 2:
                    print(f"\t\t\t- {sg} vs rest")
                    run_and_save_de(subset, sample_group_col, sg, "rest", sg_outdir)

    # ------------------------------------------------------------
    # 3) For each sample group: subset, then compare cell groups within that subset (for this cell_group_col)
    # ------------------------------------------------------------
    if has_sample_groups:
        all_sample_groups = valid_groups(adata, sample_group_col, min_cells=3)
        for sg in all_sample_groups:
            print(f"\t- Processing sample group subset for column '{cell_group_col}': {sg}")
            subset = adata[adata.obs[sample_group_col].astype(str) == sg].copy()
            ensure_categorical_str(subset, sample_group_col)
            ensure_categorical_str(subset, cell_group_col)
            cg_sample_dir = col_outdir / "sample_groups" / sanitize_filename(sg)
            groups_in_subset = valid_groups(subset, cell_group_col, min_cells=3)
            if len(groups_in_subset) < 2:
                print(f"\t\t- Skipping {cell_group_col} in sample '{sg}': fewer than 2 groups with ≥3 cells")
                continue
            for g in groups_in_subset:
                g_dir = cg_sample_dir / sanitize_filename(g)
                for other in [x for x in groups_in_subset if x != g]:
                    print(f"\t\t\t- {g} vs {other} within sample '{sg}'")
                    run_and_save_de(subset, cell_group_col, g, other, g_dir)
                if len(groups_in_subset) > 2:
                    print(f"\t\t\t- {g} vs rest within sample '{sg}'")
                    run_and_save_de(subset, cell_group_col, g, "rest", g_dir)


versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
    },
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
