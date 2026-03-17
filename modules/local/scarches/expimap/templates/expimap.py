#!/usr/bin/env python3

import os
import platform
import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scarches as sca
import pandas as pd
import scanpy as sc
import torch

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
torch.set_num_threads(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")

adata_processing = adata.copy()

if "${counts_layer}" != "X":
    adata_processing.X = adata.layers["${counts_layer}"]

# Ensure condition column exists
condition_col = "${condition_col}"
if condition_col not in adata_processing.obs.columns:
    adata_processing.obs[condition_col] = "condition"

# Prior biological knowledge in form of gene programs
if "${reference_model}":
    sca.utils.add_annotations(adata_processing, "${reference_model}", min_genes=12, clean=True)
else:
    raise ValueError("Reference model is required for EXPIMAP. Please provide a path to the reference model.")
    
# Initialization of the model with the reference network
intr_cvae = sca.models.EXPIMAP(
    adata=adata_processing,
    condition_key="${condition_col}",
    hidden_layer_sizes=${hidden_layer_sizes},
    recon_loss="${recon_loss}"
)

# Train the model
use_early_stopping = "${use_early_stopping}".lower() == 'true'
intr_cvae.train(
    n_epochs=${n_epochs},
    alpha_epoch_anneal=${alpha_epoch_anneal},
    alpha=${alpha},
    alpha_kl=${alpha_kl},
    use_early_stopping=use_early_stopping
)

# Extract the interpretable latent representation
emb = intr_cvae.get_latent(only_active=True)
adata.obsm['X_emb'] = emb  

adata.write_h5ad("${prefix}.h5ad")
df = pd.DataFrame(emb, index=adata.obs_names)
df.to_pickle("X_${prefix}.pkl")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
        "scarches": sca.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)

