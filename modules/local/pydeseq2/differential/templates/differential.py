#!/usr/bin/env python3

import os
import platform
import re

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for PyDESeq2 pseudobulk differential expression")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))


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
    counts_df = pd.DataFrame(
        np.asarray(counts),
        index=sub.obs_names,
        columns=sub.var_names,
    )

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=["donor", "condition"],
        ref_level=["condition", reference_condition],
    )
    dds.deseq2()
    safe_celltype = safe_name(celltype)

    for treatment in treatments:
        stat_res = DeseqStats(dds, contrast=["condition", str(treatment), reference_condition])
        stat_res.summary()
        out_path = f"{prefix}_{safe_celltype}_{safe_name(treatment)}_results.csv"
        stat_res.results_df.to_csv(out_path)
        written.append(out_path)

if not written:
    raise ValueError("No cell-type strata passed replicate filtering for PyDESeq2")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
