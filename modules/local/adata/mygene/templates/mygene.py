#!/usr/bin/env python3

import platform

import anndata as ad
import mygene
import yaml


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        return ad.read_zarr(path)
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(f"{stem}.zarr")


adata = read_anndata("${anndata}")

input_col = "${input_col}"

inputs = adata.var.index.to_list() if input_col == "index" else adata.var[input_col].to_list()

mg = mygene.MyGeneInfo()
mapping = {}
for i in range(0, len(inputs), 500):
    df_genes = mg.querymany(
        inputs[i : i + 500],
        scopes=["symbol", "entrezgene", "ensemblgene"],
        fields="symbol",
        species="${species}",
        as_dataframe=True,
    )
    mapping.update(df_genes["symbol"].dropna().to_dict())

outputs = [mapping.get(i, i) for i in inputs]

adata.var["symbols"] = outputs

write_anndata(adata, "${prefix}", "${anndata_format}")

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
