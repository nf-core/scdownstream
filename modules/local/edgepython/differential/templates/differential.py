#!/usr/bin/env python3

import platform
import re

import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd
import yaml

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgePython pseudobulk differential expression")


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
    counts = np.asarray(counts).T

    metadata = metadata.copy()
    metadata["condition"] = metadata["condition"].astype(str)
    metadata["donor"] = metadata["donor"].astype(str)
    metadata["condition"] = pd.Categorical(
        metadata["condition"],
        categories=[reference_condition] + treatments,
    )

    design = ep.model_matrix("~ donor + condition", metadata)
    y = ep.make_dgelist(counts=counts, samples=metadata)
    y = ep.calc_norm_factors(y)
    y = ep.estimate_disp(y)
    fit = ep.glm_ql_fit(y, design)
    safe_celltype = safe_name(celltype)

    for treatment in treatments:
        coef_name = f"condition{treatment}"
        if coef_name not in design.columns:
            continue
        coef_index = list(design.columns).index(coef_name)
        res = ep.glm_ql_ftest(fit, coef=coef_index)
        top = ep.top_tags(res, n=sub.n_vars)
        out_path = f"{prefix}_{safe_celltype}_{safe_name(treatment)}_results.csv"
        top["table"].to_csv(out_path)
        written.append(out_path)

if not written:
    raise ValueError("No cell-type strata passed replicate filtering for edgePython pseudobulk")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
