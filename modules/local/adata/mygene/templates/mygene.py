#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad
import httpx
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
try:
    df_genes = mg.querymany(
        inputs,
        scopes=["symbol", "entrezgene", "ensemblgene"],
        fields="symbol",
        species="human",
        as_dataframe=True,
    )
except httpx.HTTPStatusError as exc:
    status = exc.response.status_code
    if status >= 500:
        raise RuntimeError(
            f"mygene.info returned HTTP {status} (server error) while mapping "
            f"{len(inputs)} gene identifiers from var[{input_col!r}]. "
            "The mygene.info API is temporarily unavailable or overloaded — "
            "this is not caused by your input data. Re-run this process; "
            "if it keeps failing, check https://mygene.info or try again later."
        ) from exc
    raise RuntimeError(
        f"mygene.info returned HTTP {status} while mapping "
        f"{len(inputs)} gene identifiers from var[{input_col!r}]."
    ) from exc
except httpx.RequestError as exc:
    raise RuntimeError(
        f"Could not reach mygene.info while mapping {len(inputs)} gene identifiers "
        f"from var[{input_col!r}]: {exc}. Check network connectivity and try again."
    ) from exc
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
