#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import anndata as ad
import decoupler as dc
import pandas as pd
import yaml

adata = ad.read_h5ad("${h5ad}")
donor_col = "${donor_col}"
celltype_col = "${celltype_col}"
condition_col = "${condition_col}"
counts_layer = "${counts_layer}"
min_num_cells = int("${min_num_cells}")
min_total_counts = int("${min_total_counts}")
prefix = "${prefix}"

required_cols = [donor_col, celltype_col, condition_col]
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")

layer = None if counts_layer == "X" else counts_layer

pdata = dc.pp.pseudobulk(
    adata,
    sample_col=donor_col,
    groups_col=[celltype_col, condition_col],
    layer=layer,
    mode="sum",
)

dc.pp.filter_samples(
    pdata,
    min_cells=min_num_cells,
    min_counts=min_total_counts,
)

if pdata.n_obs == 0:
    raise ValueError("No pseudobulk samples passed filtering thresholds")

obs = pdata.obs.copy()
obs["donor"] = obs[donor_col].astype(str)
obs["celltype"] = obs[celltype_col].astype(str)
obs["condition"] = obs[condition_col].astype(str)

if "psbulk_cells" in obs.columns:
    obs["n_cells"] = obs["psbulk_cells"]
elif "psbulk_n_cells" in obs.columns:
    obs["n_cells"] = obs["psbulk_n_cells"]
else:
    raise ValueError("decoupler pseudobulk QC column psbulk_cells not found in obs")

if "psbulk_counts" in obs.columns:
    obs["total_counts"] = obs["psbulk_counts"]
else:
    raise ValueError("decoupler pseudobulk QC column psbulk_counts not found in obs")

obs["sample_id"] = obs.apply(
    lambda row: "__".join(
        str(row[col]).replace(" ", "_") for col in ("donor", "celltype", "condition")
    ),
    axis=1,
)
pdata.obs = obs
pdata.obs_names = obs["sample_id"].values

if "psbulk_props" in pdata.layers:
    del pdata.layers["psbulk_props"]

pdata.write_h5ad(f"{prefix}.h5ad")
obs[["sample_id", "donor", "celltype", "condition", "n_cells", "total_counts"]].to_csv(
    f"{prefix}_samples.tsv",
    sep="\t",
    index=False,
)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "decoupler": dc.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
