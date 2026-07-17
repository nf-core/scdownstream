#!/usr/bin/env python3

import os
import platform
import re

import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
donor_col = "${donor_col}"
condition_col = "${condition_col}"
celltype_col = "${celltype_col}"
celltype_value = "${celltype_value}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgepython_sc differential expression")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))


for col in [donor_col, condition_col, celltype_col]:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")

subset = adata[adata.obs[celltype_col].astype(str) == str(celltype_value)].copy()
if subset.n_obs < 10:
    raise ValueError(f"Too few cells ({subset.n_obs}) in cell type '{celltype_value}'")

condition_values = subset.obs[condition_col].astype(str)
conditions = sorted(condition_values.unique())
donors = sorted(subset.obs[donor_col].astype(str).unique())
if len(conditions) < 2:
    raise ValueError("At least two conditions are required for differential expression")
if len(donors) < 2:
    raise ValueError("At least two donors are required for edgepython_sc")

treatments = [condition for condition in conditions if condition != reference_condition]
if not treatments:
    raise ValueError(f"Reference condition '{reference_condition}' is the only condition present")

counts = subset.X
if sp.issparse(counts):
    counts = counts.toarray()
counts = np.asarray(counts).T

design = pd.DataFrame(
    {"Intercept": np.ones(subset.n_obs, dtype=float)},
    index=subset.obs_names,
)
for treatment in treatments:
    design[treatment] = (condition_values == treatment).astype(float)

sample_ids = subset.obs[donor_col].astype(str).to_numpy()

fit = ep.glm_sc_fit(
    counts,
    design=design,
    sample=sample_ids,
    norm_method="TMM",
)
fit = ep.shrink_sc_disp(fit, robust=True)

written = []
for treatment in treatments:
    coef = design.columns.get_loc(treatment)
    res = ep.glm_sc_test(fit, coef=coef)
    results = res["table"] if isinstance(res, dict) else res
    out_path = f"{prefix}_{safe_name(treatment)}_results.csv"
    if hasattr(results, "to_csv"):
        results.to_csv(out_path)
    else:
        pd.DataFrame(results).to_csv(out_path)
    written.append(out_path)

if not written:
    raise ValueError("No edgepython_sc contrasts could be tested")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
