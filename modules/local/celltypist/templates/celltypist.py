#!/usr/bin/env python3

import os
import platform
from pathlib import Path

os.environ["CELLTYPIST_FOLDER"] = "./tmp/celltypist"

import celltypist
import pandas as pd
import scanpy as sc
import yaml
from celltypist import models as ct_models

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

model_names = "${models.join(' ')}".split()
model_files = "${model_files}".split()

# (model name, path to a local model file or None to download it by name)
model_specs = [(name[:-4] if name.endswith(".pkl") else name, None) for name in model_names]
model_specs += [(Path(model_file).stem, model_file) for model_file in model_files]

spec_names = [name for name, _ in model_specs]
duplicate_names = sorted({name for name in spec_names if spec_names.count(name) > 1})
if duplicate_names:
    raise ValueError(f"CellTypist model names must be unique, found duplicates: {', '.join(duplicate_names)}")

adata_celltypist = adata.copy()  # make a copy of our adata
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

symbol_col = "${symbol_col}"
if symbol_col != "index" and symbol_col:
    if symbol_col not in adata_celltypist.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    adata_celltypist.var_names = adata_celltypist.var[symbol_col]

# celltypist expects a string index, because it will make unique names by appending "-1", "-2"
# to duplicates if necessary. Cast other types (e.g. CategoricalIndex) to str:
adata_celltypist.var_names = adata_celltypist.var_names.astype(str)

df_list = []
manifest_rows = []

for model_name, model_path in model_specs:
    if model_path is None:
        model_path = f"{model_name}.pkl"
        ct_models.download_models(model=model_path)
    model_obj = ct_models.Model.load(model_path)

    predictions = celltypist.annotate(adata_celltypist, model=model_obj)
    predictions_adata = predictions.to_adata()

    per_cell_col = f"annotation:celltypist:{model_name}:per_cell"
    conf_col = f"annotation:celltypist:{model_name}:per_cell:confidence"
    df_celltypist = predictions_adata.obs.loc[adata.obs.index, ["predicted_labels", "conf_score"]]
    df_celltypist.columns = [per_cell_col, conf_col]
    df_list.append(df_celltypist)
    manifest_rows.extend(
        [
            {"obs_column": per_cell_col, "aggregatable": "true"},
            {"obs_column": conf_col, "aggregatable": "false"},
        ]
    )

df_celltypist = pd.concat(df_list, axis=1)
df_celltypist.to_pickle("${prefix}.pkl")

pd.DataFrame(manifest_rows).to_csv(f"{prefix}_annotation_columns.csv", index=False)

adata.obs = pd.concat([adata.obs, df_celltypist], axis=1)
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
