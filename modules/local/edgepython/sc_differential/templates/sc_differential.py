#!/usr/bin/env python3

import os
import platform

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

for col in [donor_col, condition_col, celltype_col]:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")

subset = adata[adata.obs[celltype_col].astype(str) == str(celltype_value)].copy()
if subset.n_obs < 10:
    raise ValueError(f"Too few cells ({subset.n_obs}) in cell type '{celltype_value}'")

conditions = sorted(subset.obs[condition_col].astype(str).unique())
donors = sorted(subset.obs[donor_col].astype(str).unique())
if len(conditions) < 2:
    raise ValueError("At least two conditions are required for differential expression")
if len(donors) < 2:
    raise ValueError("At least two donors are required for edgepython_sc")

ref = reference_condition if reference_condition else conditions[0]
treatments = [condition for condition in conditions if condition != str(ref)]
if not treatments:
    raise ValueError(f"Reference condition '{ref}' is the only condition present")
treatment = treatments[0]

counts = subset.X
if sp.issparse(counts):
    counts = counts.toarray()
counts = np.asarray(counts).T

metadata = subset.obs[[donor_col, condition_col]].copy()
metadata.columns = ["donor", "condition"]
metadata["donor"] = metadata["donor"].astype(str)
metadata["condition"] = pd.Categorical(
    metadata["condition"].astype(str),
    categories=[str(ref), treatment],
)

fit = ep.glm_sc_fit(counts=counts, metadata=metadata, formula="~ condition", group_col="donor")
fit = ep.shrink_sc_disp(fit)
res = ep.glm_sc_test(fit, coef="condition", contrast=[treatment, str(ref)])
results = res["table"] if isinstance(res, dict) else res
if hasattr(results, "to_csv"):
    results.to_csv(f"{prefix}_results.csv")
else:
    pd.DataFrame(results).to_csv(f"{prefix}_results.csv")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
