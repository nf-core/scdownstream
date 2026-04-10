#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad
import scipy
import numpy as np
from scipy.sparse import csr_matrix
import platform
import yaml
import importlib.metadata

# Function borrowed from https://github.com/icbi-lab/luca/blob/5ffb0a4671e9c288b10e73de18d447ee176bef1d/lib/scanpy_helper_submodule/scanpy_helpers/util.py#L122C1-L135C21
def aggregate_duplicate_var(adata, aggr_fun=np.mean):
    retain_var = ~adata.var_names.duplicated(keep="first")
    duplicated_var = adata.var_names[adata.var_names.duplicated()].unique()
    if len(duplicated_var):
        for var in duplicated_var:
            mask = adata.var_names == var
            X_slice = adata.X[:, mask].toarray()
            var_aggr = aggr_fun(X_slice, axis=1)[:, np.newaxis]
            adata.X[:, mask] = np.repeat(var_aggr, np.sum(mask), axis=1)

        adata_dedup = adata[:, retain_var].copy()
        return adata_dedup
    else:
        return adata

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "Unknown"

    corrected = s \
        .replace(" ", "_") \
        .replace("-", "_")

    corrected = "".join([c if c.isalnum() or c == "_" else "" for c in corrected])

    # Make sure there is never more than one underscore
    corrected = corrected.replace("__", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "Unknown"

    return corrected.capitalize()

adata = ad.read_h5ad("$h5ad")

counts_layer = "${counts_layer}"
if counts_layer != "X":
    adata.X = adata.layers[counts_layer]

# Remove all obsm, varm, uns and layers
adata.obsm = {}
adata.varm = {}
adata.uns = {}
adata.layers = {}

# Convert to float32 CSR matrix
adata.X = csr_matrix(adata.X.astype(np.float32))

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
    adata.var_names = adata.var_names.str.replace(r'\\.\\d+', '', regex=True)

adata.var.index = adata.var.index.astype(str)

# Deal with duplicate genes
method = "${duplicate_var_resolution}"
if method in ["mean", "sum", "max"]:
    adata = aggregate_duplicate_var(adata, aggr_fun=getattr(np, method))
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

if label_col:
    if label_col not in adata.obs:
        raise ValueError("The specified label column does not exist in the dataset. Existing columns: " + ", ".join(adata.obs.columns))

    if label_col != "label":
        if "label" in adata.obs:
            raise ValueError("The label column already exists.")
        adata.obs["label"] = adata.obs[label_col]
        del adata.obs[label_col]

    if unknown_label != "unknown":
        if "unknown" in adata.obs["label"]:
            raise ValueError("The label column already contains 'unknown' values.")
        adata.obs["label"].replace({unknown_label: "unknown"}, inplace=True)

    # Replace all NaN values with "unknown"
    adata.obs["label"] = adata.obs["label"].astype(str)
    adata.obs["label"] = adata.obs["label"].fillna("unknown")
    adata.obs["label"] = adata.obs["label"].map(to_Florent_case)
    adata.obs["label"] = adata.obs["label"].astype("category")
else:
    if "label" in adata.obs:
        raise ValueError("The label column already exists.")
    adata.obs["label"] = "Unknown"
adata.obs["label"] = adata.obs["label"].astype("category")

# Unify conditions
condition_col = "${condition_col}"

if condition_col:
    if condition_col not in adata.obs:
        raise ValueError("The specified condition column does not exist in the dataset. Existing columns: " + ", ".join(adata.obs.columns))

    if condition_col != "condition":
        if "condition" in adata.obs:
            raise ValueError("The condition column already exists.")
        adata.obs["condition"] = adata.obs[condition_col]
        del adata.obs[condition_col]

    # Replace all NaN values with "unknown"
    adata.obs["condition"] = adata.obs["condition"].astype(str)
    adata.obs["condition"] = adata.obs["condition"].fillna("unknown")
    adata.obs["condition"] = adata.obs["condition"].map(to_Florent_case)
    adata.obs["condition"] = adata.obs["condition"].astype("category")
else:
    if "condition" in adata.obs:
        raise ValueError("The condition column already exists.")
    adata.obs["condition"] = "Unknown"
adata.obs["condition"] = adata.obs["condition"].astype("category")

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
        "anndata": importlib.metadata.version('anndata'),
        "scipy": scipy.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
