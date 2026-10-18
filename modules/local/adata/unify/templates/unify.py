#!/usr/bin/env python3

import importlib.metadata
import platform

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import yaml
from scipy.sparse import csr_matrix, hstack, issparse


def aggregate_duplicate_var(adata, method):
    groups, names = pd.factorize(adata.var_names, sort=False)
    if len(names) == adata.n_vars:
        return adata

    sizes = np.bincount(groups)
    columns = np.split(np.argsort(groups, kind="stable"), np.cumsum(sizes)[:-1])

    def aggregate(matrix):
        if not issparse(matrix):
            return np.column_stack([getattr(np, method)(matrix[:, indices], axis=1) for indices in columns])

        if method == "max":
            matrix = matrix.tocsc()
            return hstack([matrix[:, indices].max(axis=1) for indices in columns], format="csr")

        matrix = csr_matrix(matrix)
        if method == "sum" and matrix.dtype.kind in "biu":
            matrix = matrix.astype(np.uint64 if matrix.dtype.kind == "u" else np.int64, copy=False)
        weights = np.ones(adata.n_vars, dtype=matrix.dtype)
        if method == "mean":
            # The divisor includes implicit sparse zeros.
            dtype = np.float64 if matrix.dtype.kind in "biu" else matrix.dtype
            weights = weights.astype(dtype) / sizes[groups].astype(dtype)
        grouping = csr_matrix((weights, (np.arange(adata.n_vars), groups)), shape=(adata.n_vars, len(names)))
        return matrix @ grouping

    matrices = {key: aggregate(adata.layers[key]) for key in adata.layers if key is not None}
    selected = aggregate(adata.X)
    adata = adata[:, ~adata.var_names.duplicated(keep="first")].copy()
    adata.X = csr_matrix(selected, dtype=np.float32)
    for key, matrix in matrices.items():
        adata.layers[key] = matrix
    return adata


def copy_source(adata, source, role):
    matrix = adata.X if source == "X" else adata.layers.get(source)
    if matrix is None:
        raise ValueError(f"Sample {sample_id}: missing {role} counts source '{source}'")
    return matrix.copy()


def to_florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "Unknown"

    corrected = s.replace(" ", "_").replace("-", "_")

    corrected = "".join([c if c.isalnum() or c == "_" else "" for c in corrected])

    # Make sure there is never more than one underscore
    corrected = corrected.replace("__", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "Unknown"

    return corrected.capitalize()


def unify_obs_column(
    adata,
    source_col: str,
    target_name: str,
    default_value: str = "Unknown",
    *,
    apply_florent_case: bool = False,
    unknown_label: str | None = None,
) -> None:
    if source_col:
        if source_col not in adata.obs:
            raise ValueError(
                f"The specified {target_name} column does not exist in the dataset. Existing columns: "
                + ", ".join(adata.obs.columns)
            )

        if source_col != target_name:
            if target_name in adata.obs:
                raise ValueError(f"The {target_name} column already exists.")
            adata.obs[target_name] = adata.obs[source_col]
            del adata.obs[source_col]

        if unknown_label is not None and unknown_label != "unknown":
            if "unknown" in adata.obs[target_name]:
                raise ValueError(f"The {target_name} column already contains 'unknown' values.")
            adata.obs[target_name] = adata.obs[target_name].astype(str).replace({unknown_label: "unknown"})

        adata.obs[target_name] = adata.obs[target_name].astype(str)
        adata.obs[target_name] = adata.obs[target_name].fillna("unknown")
        if apply_florent_case:
            adata.obs[target_name] = adata.obs[target_name].map(to_florent_case)
    else:
        if target_name in adata.obs:
            raise ValueError(f"The {target_name} column already exists.")
        adata.obs[target_name] = default_value

    adata.obs[target_name] = adata.obs[target_name].astype("category")


adata = ad.read_h5ad("$h5ad")
sample_id = "${meta.id}"

counts_layer = "${counts_layer}"
raw_counts_layer = "${meta.raw_counts_layer ?: counts_layer}"
correction_layer = "${meta.ambient_correction_layer ?: ''}"
correction_method = "${meta.ambient_correction_method ?: 'none'}"
analysis_source = "${meta.analysis_counts_source ?: 'raw_counts'}"
canonical_correction = f"ambient_corrected_{correction_method}"

if correction_method not in {"none", "unknown", "soupx", "decontx", "cellbender", "scar"}:
    raise ValueError(f"Sample {sample_id}: invalid ambient correction method '{correction_method}'")
if bool(correction_layer) != (correction_method not in {"none", "unknown"}):
    raise ValueError(
        f"Sample {sample_id}: correction layer '{correction_layer}' is inconsistent with method '{correction_method}'"
    )
if analysis_source != "raw_counts" and not (correction_layer and analysis_source == canonical_correction):
    raise ValueError(
        f"Sample {sample_id}: analysis counts source '{analysis_source}' has no generated correction layer"
    )
expected_source = raw_counts_layer if analysis_source == "raw_counts" else correction_layer
if counts_layer != expected_source:
    raise ValueError(
        f"Sample {sample_id}: selected counts layer '{counts_layer}' does not match analysis source '{expected_source}'"
    )

# Reserved destinations can alias any input source.
selected = copy_source(adata, counts_layer, "selected analysis")
retained = {"raw_counts": copy_source(adata, raw_counts_layer, "raw")}
if correction_layer:
    retained[canonical_correction] = copy_source(adata, correction_layer, "ambient correction")

adata.obsm = {}
adata.varm = {}
adata.uns = {}
adata.layers = retained

adata.X = csr_matrix(selected, dtype=np.float32)
adata.obs["ambient_correction_method"] = pd.Categorical([correction_method] * adata.n_obs)
adata.obs["analysis_counts_source"] = pd.Categorical([analysis_source] * adata.n_obs)

# Unify gene symbols
symbol_col = "${symbol_col}"

label_col = "${label_col}"
unknown_label = "${unknown_label}"

if symbol_col != "index":
    assert symbol_col in adata.var.columns, f"Symbol column {symbol_col} not found in var table"
    adata.var["original_index"] = adata.var.index
    adata.var.index = adata.var[symbol_col]

if "${aggregate_isoforms}" == "true":
    # Remove all numeric suffixes following a dot, keep non-numeric suffixes
    adata.var_names = adata.var_names.str.replace(r"\\.\\d+", "", regex=True)

adata.var.index = adata.var.index.astype(str)

# Deal with duplicate genes
method = "${duplicate_var_resolution}"
if method in ["mean", "sum", "max"]:
    adata = aggregate_duplicate_var(adata, method)
elif method == "make_unique":
    adata.var_names_make_unique()
else:
    raise ValueError(f"Invalid aggregation method: {method}")

# Prevent duplicate cells
adata.obs_names_make_unique()
adata.obs_names = "${meta.id}_" + adata.obs_names

# Unify batches
batch_col = "${batch_col}"
if batch_col not in adata.obs:
    adata.obs[batch_col] = "${meta.id}"

if batch_col != "batch":
    if "batch" in adata.obs:
        raise ValueError("The batch column already exists.")
    adata.obs["batch"] = adata.obs[batch_col]
    if batch_col not in [label_col, "sample"]:
        del adata.obs[batch_col]
adata.obs["batch"] = adata.obs["batch"].astype(str).astype("category")

unify_obs_column(
    adata,
    label_col,
    "label",
    default_value="Unknown",
    apply_florent_case=True,
    unknown_label=unknown_label,
)

# Unify conditions
condition_col = "${condition_col}"

unify_obs_column(
    adata,
    condition_col,
    "condition",
    default_value="Unknown",
    apply_florent_case=True,
)

# Unify donor / biological replicate column for pseudobulk analyses
donor_col = "${donor_col}"

unify_obs_column(
    adata,
    donor_col,
    "donor",
    default_value="unknown",
)

# Add "sample" column
if "sample" in adata.obs and not adata.obs["sample"].equals("${meta.id}"):
    adata.obs["sample_original"] = adata.obs["sample"]
adata.obs["sample"] = "${meta.id}"
adata.obs["sample"] = adata.obs["sample"].astype("category")

# Add sample to batch column, to avoid overlap with other samples
adata.obs["batch"] = adata.obs["batch"].astype(str) + "_" + adata.obs["sample"].astype(str)
adata.obs["batch"] = adata.obs["batch"].astype("category")

# Ensure no named indices remain (nft-anndata doesn't support them)
# TODO: Remove this after nft-anndata plugin is updated to support named indices
if adata.var.index.name is not None:
    adata.var.index.name = None
if adata.obs.index.name is not None:
    adata.obs.index.name = None

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "scipy": scipy.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
