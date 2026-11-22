#!/usr/bin/env python3

import platform
import re
from pathlib import Path

import anndata as ad
import edgepython as ep
import h5py
import numpy as np
import pandas as pd
import yaml

prefix = "${prefix}"
donor_col = "${donor_col}"
condition_col = "${condition_col}"
celltype_col = "${celltype_col}"
celltype_value = "${celltype_value}"
reference_condition = "${reference_condition}"
if not reference_condition:
    raise ValueError("reference_condition must be set for edgepython_sc differential expression")

STANDARD_DE_COLUMNS = ["gene", "log2fc", "pvalue", "padj", "group", "contrast", "stratum"]


def write_standard_de_parquet(df, path, *, genes, group, contrast, stratum):
    out = df.rename(columns={"logFC": "log2fc", "PValue": "pvalue", "FDR": "padj"})
    out["gene"] = np.asarray(genes, dtype=str)
    for col in ("log2fc", "pvalue", "padj"):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out["group"] = str(group)
    out["contrast"] = str(contrast)
    out["stratum"] = str(stratum)
    extra = [col for col in out.columns if col not in STANDARD_DE_COLUMNS]
    # Absolute path so pyarrow does not treat colons in Nextflow prefixes as URI schemes.
    out[STANDARD_DE_COLUMNS + extra].to_parquet(str(Path(path).resolve()), index=False, engine="pyarrow")


with h5py.File("${h5ad}", "r") as handle:
    obs = ad.io.read_elem(handle["obs"])
    for col in (donor_col, condition_col, celltype_col):
        if col not in obs.columns:
            raise ValueError(f"Column '{col}' not found in adata.obs")

    rows = np.flatnonzero(obs[celltype_col].astype(str).to_numpy() == str(celltype_value))
    if len(rows) < 10:
        raise ValueError(f"Too few cells ({len(rows)}) in cell type '{celltype_value}'")

    # edgePython groups donors in first-seen order; pre-grouping avoids its dense copy.
    sample_ids = obs[donor_col].iloc[rows].astype(str).to_numpy()
    donor_codes, unique_donors = pd.factorize(sample_ids, sort=False)
    if len(unique_donors) < 2:
        raise ValueError("At least two donors are required for edgepython_sc")

    order = np.argsort(donor_codes, kind="stable")
    rows = rows[order]
    sample_ids = sample_ids[order]
    condition_values = obs[condition_col].iloc[rows].astype(str).to_numpy()

    conditions = sorted(set(condition_values))
    if len(conditions) < 2:
        raise ValueError("At least two conditions are required for differential expression")
    treatments = [c for c in conditions if c != reference_condition]

    gene_names = np.asarray(ad.io.read_elem(handle["var"][handle["var"].attrs["_index"]]), dtype=str)

    matrix = handle["X"]
    if isinstance(matrix, h5py.Group):
        counts = ad.io.sparse_dataset(matrix)[rows, :].astype(np.float64).toarray().T
    else:
        counts = np.empty((len(rows), matrix.shape[1]), dtype=np.float64)
        for start in range(0, len(rows), 256):
            # HDF5 dense indexing requires increasing row indices.
            indices = start + np.argsort(rows[start : start + 256])
            counts[indices] = matrix[rows[indices], :]
        counts = counts.T

design = np.column_stack([np.ones(len(rows)), *(condition_values == t for t in treatments)])

fit = ep.glm_sc_fit(counts, design=design, sample=sample_ids, norm_method="TMM")
fit = ep.shrink_sc_disp(fit, robust=True)

gene_mask = fit.get("gene_mask")
gene_labels = gene_names[np.asarray(gene_mask, dtype=bool)] if gene_mask is not None else gene_names

for coef, treatment in enumerate(treatments, start=1):
    results = ep.glm_sc_test(fit, coef=coef)["table"]
    # glm_sc_test returns coefficients on the natural-log link scale.
    results["logFC"] = results["logFC"] / np.log(2)
    safe_treatment = re.sub(r"[^A-Za-z0-9._-]+", "_", treatment)
    write_standard_de_parquet(
        results,
        f"{prefix}_{safe_treatment}_results.parquet",
        genes=gene_labels,
        group=treatment,
        contrast=f"{treatment} vs {reference_condition}",
        stratum=f"celltype={celltype_value}",
    )

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "edgepython": ep.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
