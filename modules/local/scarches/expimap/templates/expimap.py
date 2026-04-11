#!/usr/bin/env python3

import argparse
import os
import shlex
import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad

# Monkey-patch anndata for scarches compatibility
ad.read = ad.read_h5ad

import scarches as sca
import pandas as pd
import scanpy as sc
import torch

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
torch.set_num_threads(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
args_str = "${args}"

parser = argparse.ArgumentParser()
parser.add_argument(
    "--hidden-layer-sizes",
    nargs="+",
    type=int,
    default=[256, 256, 256],
)
parser.add_argument(
    "--recon-loss",
    default="nb",
    choices=["nb", "zinb", "poisson", "normal"],
)
parser.add_argument("--n-epochs", type=int, default=400)
parser.add_argument("--alpha-epoch-anneal", type=int, default=100)
parser.add_argument("--alpha", type=float, default=0.7)
parser.add_argument("--alpha-kl", type=float, default=0.5)
parser.add_argument(
    "--no-use-early-stopping",
    dest="use_early_stopping",
    action="store_false",
)
parser.set_defaults(use_early_stopping=True)
cli = parser.parse_args(shlex.split(args_str) if args_str.strip() else [])

adata_processing = adata.copy()

if "${counts_layer}" != "X":
    adata_processing.X = adata.layers["${counts_layer}"]

raw_condition = "${condition_col}".strip()
if not raw_condition:
    raise ValueError(
        "EXPIMAP requires a non-empty condition_col: name an existing adata.obs column "
        "(e.g. batch or condition). scArches EXPIMAP does not support training without a condition key."
    )
if raw_condition not in adata_processing.obs.columns:
    raise ValueError(
        f"EXPIMAP condition_col {raw_condition!r} is not present in adata.obs.columns"
    )
condition_key = raw_condition

if "${reference_model}":
    sca.utils.add_annotations(
        adata_processing, "${reference_model}", min_genes=12, clean=True
    )
else:
    raise ValueError(
        "Reference model is required for EXPIMAP. Please provide a path to the reference model."
    )

intr_cvae = sca.models.EXPIMAP(
    adata=adata_processing,
    condition_key=condition_key,
    hidden_layer_sizes=list(cli.hidden_layer_sizes),
    recon_loss=cli.recon_loss,
)

intr_cvae.train(
    n_epochs=cli.n_epochs,
    alpha_epoch_anneal=cli.alpha_epoch_anneal,
    alpha=cli.alpha,
    alpha_kl=cli.alpha_kl,
    use_early_stopping=cli.use_early_stopping,
)

emb = intr_cvae.get_latent(only_active=True)
adata.obsm["X_emb"] = emb

adata.write_h5ad(f"{prefix}.h5ad")
df = pd.DataFrame(emb, index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

versions = {
    "${task.process}": {
        "scarches": sca.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
