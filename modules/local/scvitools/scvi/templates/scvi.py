#!/usr/bin/env python3

import scvi
import anndata as ad
import pandas as pd
from scvi.model import SCVI
import platform
import torch

torch.set_float32_matmul_precision('medium')

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
scvi.settings.num_threads = int("${task.cpus}")

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

adata = ad.read_h5ad("${h5ad}")

setup_kwargs = {"batch_key": "batch"}

# Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
model_kwargs = {
    "gene_likelihood": "nb",
    "n_layers": 2,
    "n_hidden": 128,
    "n_latent": 30,
}

train_kwargs = {"train_size": 1.0}

n_epochs = int(min([round((20000 / adata.n_obs) * 400), 400]))

SCVI.setup_anndata(adata, **setup_kwargs)
model = SCVI(adata, **model_kwargs)

if "${task.ext.use_gpu}" == "true":
    model.to_device(0)

model.train(max_epochs=n_epochs, **train_kwargs)

adata.obsm["X_emb"] = model.get_latent_representation()

adata.write_h5ad("${prefix}.h5ad")
model.save("${prefix}_model")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle("X_${prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scvi": scvi.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
