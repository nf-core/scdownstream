#!/usr/bin/env python3

import os
from pathlib import Path

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ.setdefault("TORCHINDUCTOR_CACHE_DIR", os.path.join(os.getcwd(), "torch_cache"))

import anndata as ad
import pyarrow
import scvi
import torch
import yaml
from scvi.external import SOLO
from scvi.model import SCVI
from threadpoolctl import threadpool_limits

torch.set_float32_matmul_precision("medium")
scvi.settings.seed = 0

threadpool_limits(int("${task.cpus}"))
scvi.settings.num_threads = int("${task.cpus}")

batch_key = "${batch_key ?: ''}"
max_epochs = "${max_epochs ?: ''}"


def train_model(model):
    if "${task.accelerator ? 'true' : 'false'}" == "true":
        model.to_device(0)

    model.train(max_epochs=int(max_epochs) if max_epochs else None, datasplitter_kwargs={"drop_last": True})


adata = ad.read_h5ad("${h5ad}")

if batch_key:
    SCVI.setup_anndata(adata, batch_key=batch_key)
else:
    SCVI.setup_anndata(adata)

scvi_model = SCVI(adata)

train_model(scvi_model)

adata.obs["doublet"] = False

batches = adata.obs[batch_key].unique() if batch_key else [0]
for batch in batches:
    solo_model = SOLO.from_scvi_model(scvi_model, restrict_to_batch=batch if len(batches) > 1 else None)

    train_model(solo_model)
    result = solo_model.predict(False)

    doublets = result[result == "doublet"].index.tolist()
    adata.obs.loc[doublets, "doublet"] = True


df = adata.obs[["doublet"]].copy()
df.columns = ["${prefix}"]
if df.index.hasnans or not df.index.is_unique or not all(isinstance(cell, str) for cell in df.index):
    raise ValueError("Doublet predictions require unique, non-null string cell IDs")
df.insert(0, "__index__", df.index)
df.to_parquet(Path("${prefix}.parquet").resolve(), index=False, engine="pyarrow", compression="snappy")

adata = adata[~adata.obs["doublet"]].copy()
adata.obs.drop("doublet", axis=1, inplace=True)

del adata.uns["_scvi_manager_uuid"]
del adata.uns["_scvi_uuid"]

adata.write_h5ad("${prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "scvi": scvi.__version__,
        "pyarrow": pyarrow.__version__,
        "torch": str(torch.__version__),
        "cuda": torch.version.cuda or "no CUDA available",
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))
