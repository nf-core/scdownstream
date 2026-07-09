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

    ref = reference_condition if reference_condition else conditions[0]
    treatments = [condition for condition in conditions if condition != str(ref)]
    if not treatments:
        continue
    treatment = treatments[0]

    counts = sub.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    counts_df = pd.DataFrame(np.asarray(counts).T, index=sub.var_names, columns=sub.obs_names)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=["donor", "condition"],
        ref_level=["condition", str(ref)],
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", str(treatment), str(ref)])
    stat_res.summary()
    safe_celltype = re.sub(r"[^A-Za-z0-9._-]+", "_", str(celltype))
    out_path = f"{prefix}_{safe_celltype}_results.csv"
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
