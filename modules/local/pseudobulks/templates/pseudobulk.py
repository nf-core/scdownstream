#!/usr/bin/env python3

import os
import platform
import random
import sys

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import pandas as pd
import numpy as np
import scipy

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

#read and harmonize parameters:
sample_key = "batch" #pipeline's sample key, default pseudobulking column
cell_identity_key="label" #pipeline's cell_type key

donors_to_drop = []

adata = sc.read_h5ad("${h5ad}")

if "${groups}" != "null": # column names to consider for pseudobulking
    group_keys = "${groups}".replace(' ','').split(',') # 1-n columns
else:
    group_keys = [sample_key] # no list given

mode = "${mode}" #sum, mean, median

if "${min_cells}" != "null": # minimal cells per donor
    min_cells = int("${min_cells}")
else:
    min_cells = 10

if "${replicates_per_sample}" != "null":
    replicates_per_patient = int("${replicates_per_sample}")
    pseudoreplicates = True
else:
    replicates_per_patient = 1
    pseudoreplicates = False

#create bulk column:
if len(group_keys)>1:
    adata.obs['bulk_column'] = adata.obs[group_keys].agg('_'.join, axis=1)
    sample_key = 'bulk_column' # overwrite sample_key

elif group_keys is not None: # needed to specifiy to user: in case no groups are given, pseudobulks will be built on samples/batch
    sample_key = group_keys[0] # overwrite sample key 

# best practice: 30 cells per donor minimum here 10 enough:
size_by_donor = adata.obs.groupby([sample_key]).size()
donors_to_drop = [
    donor for donor in size_by_donor.index
    if size_by_donor[donor] < min_cells
]

has_enough_donors = (size_by_donor >= min_cells).any()

if not has_enough_donors:
    sys.exit(f"Error: No bulks have more than {min_cells} cells. Please provide less bulking groups or more samples.")

adata.obs[sample_key] = adata.obs[sample_key].astype("category") # needed for pseudoreplicates

df = pd.DataFrame(columns=[*adata.var_names,])


if pseudoreplicates:
    for i, donor in enumerate(donors := adata.obs[sample_key].cat.categories): # from best practices book
        if donor not in donors_to_drop:
            adata_donor = adata[adata.obs[sample_key] == donor]
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: mode for gene in adata_replicate.var_names}
                # create a df with all genes, donor and condition info
                # Convert X to dense if it's sparse
                X_dense = adata_replicate.X.toarray() if scipy.sparse.issparse(adata_replicate.X) else adata_replicate.X
                df_donor = pd.DataFrame(X_dense)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                df_donor = df_donor.join(adata_replicate.obs)
                # aggregate
                df_donor = df_donor.groupby(sample_key).agg(agg_dict)
                df_donor[sample_key] = donor
                df.loc[f"{donor}_{i}"] = df_donor.loc[donor]
                adata[adata.obs[sample_key] == donor].obs['psbulk'] = df_donor.loc[donor]
    adata.varm['psbulk'] = df.T.astype(str)

else:
    for sample in adata.obs[sample_key].unique(): # from best practices book and https://www.youtube.com/watch?v=Ee0PQUwVH8Q, skipping pseudoreplicates
        if sample in donors_to_drop:
            continue
        sample_subset = adata[adata.obs[sample_key] == sample]
        agg_dict = {gene: mode for gene in sample_subset.var_names}
        X_dense = sample_subset.X.toarray() if scipy.sparse.issparse(sample_subset.X) else sample_subset.X
        df_donor = pd.DataFrame(X_dense)
        df_donor.index = sample_subset.obs_names
        df_donor.columns = sample_subset.var_names
        df_donor = df_donor.join(sample_subset.obs)
        df_donor = df_donor.groupby(sample_key).agg(agg_dict)
        df_donor[sample_key] = sample
        df.loc[f"{sample}"] = df_donor.loc[sample]
        adata[adata.obs[sample_key] == sample].obs['psbulk'] = df_donor.loc[sample]
        
    adata.varm['psbulk'] = df.T.astype(str)

# save df
df.to_pickle("${prefix}.pkl")

#save adata
adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
