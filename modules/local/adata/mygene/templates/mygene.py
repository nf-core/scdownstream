#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad
import mygene
import yaml

adata = ad.read_h5ad("$h5ad")

input_col = "${input_col}"
output_col = "${output_col}"

inputs = (
    adata.var.index.to_list()
    if input_col == "index"
    else adata.var[input_col].to_list()
)

mg = mygene.MyGeneInfo()
df_genes = mg.querymany(inputs,
    scopes=["symbol", "entrezgene", "ensemblgene"],
    fields="symbol", species="human", as_dataframe=True)
mapping = df_genes["symbol"].dropna().to_dict()

outputs = [mapping.get(i, i) for i in inputs]

if output_col == "index":
    adata.var.index = outputs
else:
    adata.var[output_col] = outputs

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "mygene": mygene.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
