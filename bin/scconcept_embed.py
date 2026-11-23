#!/usr/bin/env python3

import argparse
import hashlib
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from concept import scConcept
from scipy.sparse import issparse

parser = argparse.ArgumentParser(description="Embed human cells with the pinned scConcept checkpoint.")
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--model", type=Path, required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--threads", type=int, required=True)
args = parser.parse_args()

bundle = args.model
hashes = {
    "config.yaml": "fafabeccbd5693e7dcc48097e3ad8175aecab849d440f805be71b669ac86497f",
    "model.ckpt": "ed4a51d77b3205029ca7253d225768761780ef620aad98ed496c62ba59f3fe29",
    "gene_mappings/hsapiens.csv": "fed5b8da96c277d09049b8922a6aaa432deee75a159cd6c14128129110d5ef8c",
}
for relative_path, expected in hashes.items():
    with (bundle / relative_path).open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != expected:
            raise ValueError(
                f"scConcept model file {relative_path} does not match the pinned corpus40M-model30M release."
            )

adata = ad.read_h5ad(args.input)
if adata.n_obs == 0:
    raise ValueError("scConcept requires at least one cell.")

vocabulary = pd.read_csv(bundle / "gene_mappings/hsapiens.csv", dtype=str)
vocabulary = vocabulary[vocabulary["gene_id"].str.startswith("ENSG")]
vocabulary_ids = set(vocabulary["gene_id"])
symbols = vocabulary.dropna(subset=["gene_name"]).assign(symbol=lambda df: df["gene_name"].str.casefold())
symbol_to_id = symbols.drop_duplicates("symbol", keep=False).set_index("symbol")["gene_id"].to_dict()

ensg = re.compile(r"(ENSG[0-9]+)(?:[.][0-9]+)?")
gene_ids = pd.Index(
    [match[1] if (match := ensg.fullmatch(name)) else symbol_to_id.get(name.casefold()) for name in adata.var_names]
)
matched = gene_ids.isin(vocabulary_ids)
print(f"scConcept vocabulary coverage: {matched.sum()} of {adata.n_vars} input genes matched.", flush=True)
if not matched.any():
    raise ValueError("No genes matched the human scConcept vocabulary. Supply human gene symbols or Ensembl gene IDs.")
if gene_ids[matched].has_duplicates:
    raise ValueError(
        "Multiple input genes map to the same scConcept Ensembl ID. Resolve duplicate genes before embedding."
    )

counts = adata[:, matched].X
values = counts.data if issparse(counts) else counts
if not np.isfinite(values).all() or (values < 0).any():
    raise ValueError("scConcept requires finite, non-negative counts in AnnData X.")
empty = np.asarray((counts > 0).sum(axis=1)).ravel() == 0
if empty.any():
    raise ValueError(f"{empty.sum()} cells have no expression among genes in the scConcept vocabulary.")

model_input = ad.AnnData(X=counts, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=gene_ids[matched]))
torch.set_num_threads(args.threads)
concept = scConcept()
concept.load_config_and_model(
    config=bundle / "config.yaml",
    model_path=bundle / "model.ckpt",
    gene_mappings_path=bundle / "gene_mappings",
)
output = concept.extract_embeddings(
    model_input,
    species="hsapiens",
    batch_size=16,
    num_workers=0,
    return_type="torch",
)
embedding = output["cls_cell_emb"].float().cpu().numpy()
if embedding.shape != (adata.n_obs, 512) or not np.isfinite(embedding).all():
    raise ValueError("scConcept did not produce a finite 512-dimensional embedding for each input cell.")

adata.obsm["X_emb"] = embedding
adata.write_h5ad(f"{args.prefix}.h5ad", compression="gzip")
table = pd.DataFrame(embedding, index=adata.obs_names, columns=[str(i) for i in range(embedding.shape[1])])
table.reset_index(names="__index__").to_parquet(f"X_{args.prefix}.parquet", engine="pyarrow", index=False)
