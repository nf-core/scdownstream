#!/usr/bin/env python3

import argparse
import shlex
from pathlib import Path

import anndata as ad
import yaml

# Monkey-patch anndata for scarches compatibility
ad.read = ad.read_h5ad

import pandas as pd
import pyarrow as pa
import scanpy as sc
import scarches as sca
import torch
from threadpoolctl import threadpool_limits


def write_fragment(frame, path):
    fragment = frame.copy()
    fragment.columns = fragment.columns.map(str)
    if "__index__" in fragment.columns or not fragment.columns.is_unique:
        raise ValueError("Fragment columns must be unique and must not contain reserved column '__index__'.")
    if (
        fragment.index.nlevels != 1
        or not fragment.index.is_unique
        or fragment.index.hasnans
        or not all(isinstance(value, str) for value in fragment.index)
    ):
        raise ValueError("Fragment axis identifiers must be unique, non-null strings.")
    fragment = fragment.reset_index(names="__index__")
    fragment.to_parquet(str(Path(path).resolve()), engine="pyarrow", index=False, compression="snappy")


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
    default=None,
)
parser.add_argument(
    "--recon-loss",
    default=None,
    choices=["mse", "nb"],
)
parser.add_argument("--n-epochs", type=int, default=None)
parser.add_argument("--alpha-epoch-anneal", type=int, default=None)
parser.add_argument("--alpha", type=float, default=None)
parser.add_argument("--alpha-kl", type=float, default=None)
early = parser.add_mutually_exclusive_group()
early.add_argument(
    "--use-early-stopping",
    dest="use_early_stopping",
    action="store_const",
    const=True,
    default=None,
)
early.add_argument(
    "--no-use-early-stopping",
    dest="use_early_stopping",
    action="store_const",
    const=False,
    default=None,
)
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
    raise ValueError(f"EXPIMAP condition_col {raw_condition!r} is not present in adata.obs.columns")
condition_key = raw_condition

if "${reference_model}":
    sca.utils.add_annotations(adata_processing, "${reference_model}", min_genes=12, clean=True)
else:
    raise ValueError("Reference model is required for EXPIMAP. Please provide a path to the reference model.")

expimap_kw = {
    "adata": adata_processing,
    "condition_key": condition_key,
}
if cli.hidden_layer_sizes is not None:
    expimap_kw["hidden_layer_sizes"] = list(cli.hidden_layer_sizes)
if cli.recon_loss is not None:
    expimap_kw["recon_loss"] = cli.recon_loss
intr_cvae = sca.models.EXPIMAP(**expimap_kw)

train_kw = {}
for name in (
    "n_epochs",
    "alpha_epoch_anneal",
    "alpha",
    "alpha_kl",
    "use_early_stopping",
):
    val = getattr(cli, name)
    if val is not None:
        train_kw[name] = val
intr_cvae.train(**train_kw)

emb = intr_cvae.get_latent(only_active=True)
adata.obsm["X_emb"] = emb

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
df = pd.DataFrame(emb, index=adata.obs_names)
write_fragment(df, f"X_{prefix}.parquet")

versions = {
    "${task.process}": {
        "scarches": sca.__version__,
        "pyarrow": pa.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
