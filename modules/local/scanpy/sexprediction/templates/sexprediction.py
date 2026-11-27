#!/usr/bin/env python3

import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
import scanpy as sc
import yaml
from scipy import sparse
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

symbol_col = "${symbol_col}"
prefix = "${prefix}"
sample_id = "${meta.id}"

# A sample is called male or female when at least this share of its informative
# cells agree. Informative cells express Y genes or X-inactivation genes, not both.
CALL_THRESHOLD = 0.8
# Samples with fewer informative cells than this are reported as unknown.
MIN_INFORMATIVE_CELLS = 20


def read_gene_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def summed_counts(matrix):
    totals = matrix.sum(axis=1)
    return np.asarray(totals, dtype=np.float64).ravel()


y_genes = read_gene_list("${y_genes}")
x_genes = read_gene_list("${x_genes}")

adata = sc.read_h5ad("${h5ad}")

if symbol_col == "index":
    symbols = pd.Index(adata.var_names.astype(str))
else:
    if symbol_col not in adata.var.columns:
        raise ValueError(
            f"symbol_col '{symbol_col}' not found in adata.var. Available columns: {list(adata.var.columns)}"
        )
    symbols = pd.Index(adata.var[symbol_col].astype(str))

y_mask = symbols.isin(y_genes)
x_mask = symbols.isin(x_genes)
if not y_mask.any():
    print(f"Warning: none of the Y genes {y_genes} were found in {sample_id}; male cells cannot be detected.")
if not x_mask.any():
    print(
        f"Warning: none of the X-inactivation genes {x_genes} were found in {sample_id}; female cells cannot be detected."
    )

matrix = adata.X
if sparse.issparse(matrix):
    matrix = matrix.tocsr()
    values = matrix.data
else:
    values = np.asarray(matrix)
if values.size and not np.allclose(values, np.round(values)):
    print(
        f"Warning: adata.X of {sample_id} does not look like raw counts; fractions are computed on the values as given."
    )

total_counts = summed_counts(matrix)
y_counts = summed_counts(matrix[:, np.flatnonzero(y_mask)]) if y_mask.any() else np.zeros(adata.n_obs)
x_counts = summed_counts(matrix[:, np.flatnonzero(x_mask)]) if x_mask.any() else np.zeros(adata.n_obs)

with np.errstate(divide="ignore", invalid="ignore"):
    y_fraction = np.where(total_counts > 0, y_counts / total_counts, 0.0)
    x_fraction = np.where(total_counts > 0, x_counts / total_counts, 0.0)

has_y = y_counts > 0
has_x = x_counts > 0
cell_call = np.select(
    [has_y & ~has_x, has_x & ~has_y, has_x & has_y],
    ["male", "female", "both"],
    default="unknown",
)

n_male = int((cell_call == "male").sum())
n_female = int((cell_call == "female").sum())
n_informative = n_male + n_female
male_share = n_male / n_informative if n_informative else float("nan")

if n_informative < MIN_INFORMATIVE_CELLS:
    predicted_sex = "unknown"
elif male_share >= CALL_THRESHOLD:
    predicted_sex = "male"
elif male_share <= 1 - CALL_THRESHOLD:
    predicted_sex = "female"
else:
    predicted_sex = "mixed"

adata.obs["sex_y_fraction"] = y_fraction
adata.obs["sex_xist_fraction"] = x_fraction
adata.obs["sex_cell_call"] = pd.Series(cell_call, index=adata.obs_names, dtype=str)
adata.obs["predicted_sex"] = pd.Series(predicted_sex, index=adata.obs_names, dtype=str)

columns = ["sex_y_fraction", "sex_xist_fraction", "sex_cell_call", "predicted_sex"]
fragment = adata.obs[columns].copy()
fragment["sex_cell_call"] = fragment["sex_cell_call"].astype(str)
fragment["predicted_sex"] = fragment["predicted_sex"].astype(str)
fragment.reset_index(names="__index__").to_parquet(
    Path(f"{prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
)
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

total = total_counts.sum()
summary = pd.DataFrame(
    [
        {
            "sample": sample_id,
            "n_cells": adata.n_obs,
            "n_male_cells": n_male,
            "n_female_cells": n_female,
            "n_both_cells": int((cell_call == "both").sum()),
            "n_unknown_cells": int((cell_call == "unknown").sum()),
            "male_share": male_share,
            "y_fraction": y_counts.sum() / total if total else 0.0,
            "xist_fraction": x_counts.sum() / total if total else 0.0,
            "predicted_sex": predicted_sex,
        }
    ]
)
summary.to_csv(f"{prefix}.tsv", sep="\\t", index=False)

print(f"{sample_id}: predicted sex {predicted_sex} ({n_male} male and {n_female} female informative cells)")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
