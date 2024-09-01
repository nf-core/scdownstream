#!/usr/bin/env python3

import scanpy as sc
import scipy
import numpy as np
from scipy.sparse import csr_matrix
import platform

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

# Function borrowed from https://github.com/icbi-lab/luca/blob/5ffb0a4671e9c288b10e73de18d447ee176bef1d/lib/scanpy_helper_submodule/scanpy_helpers/util.py#L122C1-L135C21
def aggregate_duplicate_var(adata, aggr_fun=np.mean):
    retain_var = ~adata.var_names.duplicated(keep="first")
    duplicated_var = adata.var_names[adata.var_names.duplicated()].unique()
    if len(duplicated_var):
        for var in duplicated_var:
            mask = adata.var_names == var
            var_aggr = aggr_fun(adata.X[:, mask], axis=1)[:, np.newaxis]
            adata.X[:, mask] = np.repeat(var_aggr, np.sum(mask), axis=1)

        adata_dedup = adata[:, retain_var].copy()
        return adata_dedup
    else:
        return adata

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "unknown"

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
        return "unknown"

    return corrected.capitalize()

adata = sc.read_h5ad("$h5ad")

# Convert to float32 CSR matrix
adata.X = csr_matrix(adata.X.astype(np.float32))

# Prevent duplicate cells
adata.obs_names_make_unique()
adata.obs_names = "${meta.id}_" + adata.obs_names

# Remove all obsm, varm, uns and layers
adata.obsm = {}
adata.varm = {}
adata.layers = {}
adata.uns = {}

# Unify batches
batch_col = "${meta.batch_col}"
if batch_col not in adata.obs:
    adata.obs[batch_col] = "${meta.id}"

if batch_col != "batch":
    if "batch" in adata.obs:
        raise ValueError("The batch column already exists.")
    adata.obs["batch"] = adata.obs[batch_col]
    del adata.obs[batch_col]

# Unify labels
label_col = "${meta.label_col ?: ''}"
unknown_label = "${meta.unknown_label}"

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
    adata.obs["label"] = "unknown"

# Unify gene symbols
symbol_col = "${meta.symbol_col ?: 'index'}"

if symbol_col != "index":
    if symbol_col == "none":
        import mygene
        mg = mygene.MyGeneInfo()
        df_genes = mg.querymany(adata.var.index,
            scopes=["symbol", "entrezgene", "ensemblgene"],
            fields="symbol", species="human", as_dataframe=True)
        mapping = df_genes["symbol"].dropna().to_dict()

        adata.var.index = adata.var.index.map(lambda x: mapping.get(x, x))
    else:
        adata.var.index = adata.var[symbol_col]
        del adata.var[symbol_col]

# Aggregate duplicate genes
method = "${params.var_aggr_method}"
if not method in ["mean", "sum", "max"]:
    raise ValueError(f"Invalid aggregation method: {method}")

adata = aggregate_duplicate_var(adata, aggr_fun=getattr(np, method))

# Add "sample" column
if "sample" in adata.obs and not adata.obs["sample"].equals("${meta.id}"):
    adata.obs["sample_original"] = adata.obs["sample"]
adata.obs["sample"] = "${meta.id}"

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "scipy": scipy.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
