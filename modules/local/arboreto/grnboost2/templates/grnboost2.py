#!/usr/bin/env python3

import multiprocessing
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
import scanpy as sc
import yaml
from arboreto import core
from scipy import sparse
from threadpoolctl import threadpool_limits

cpus = int("${task.cpus}")
prefix = "${prefix}"
obs_key = "${obs_key}"
tfs_path = "${tfs}"
n_top_genes = int("${n_top_genes}")
min_cells = int("${min_cells}")
max_links = int("${max_links}")
seed = int("${seed}")

COLUMNS = ["group", "TF", "target", "importance"]


def read_gene_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def dense(matrix):
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


adata = sc.read_h5ad("${h5ad}")
adata.var_names = adata.var_names.astype(str)
if obs_key not in adata.obs.columns:
    raise ValueError(f"obs_key '{obs_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")

# Target genes are the most variable genes across all cells, selected on log-normalised data.
# The network itself is inferred from the counts in X, as in GRNBoost2 and pySCENIC.
lognorm = adata.copy()
sc.pp.normalize_total(lognorm)
sc.pp.log1p(lognorm)
sc.pp.highly_variable_genes(lognorm, n_top_genes=min(n_top_genes, adata.n_vars))
target_genes = adata.var_names[lognorm.var["highly_variable"].to_numpy()]
del lognorm

if tfs_path:
    tf_list = read_gene_list(tfs_path)
    regulators = adata.var_names[adata.var_names.isin(tf_list)]
    if len(regulators) == 0:
        raise ValueError(f"None of the {len(tf_list)} genes in {tfs_path} were found in adata.var_names")
    print(f"Using {len(regulators)} of {len(tf_list)} listed regulators found in the data.")
else:
    regulators = target_genes
    print("No regulator list given; using all selected target genes as candidate regulators.")

genes = target_genes.union(regulators, sort=False)
target_indices = genes.get_indexer(target_genes)
print(f"Inferring networks for {len(target_genes)} target genes with {len(regulators)} candidate regulators.")

# State shared with worker processes through fork
_tf_matrix = None
_tf_names = None
_expression = None


def infer_target(index):
    target = genes[index]
    target_expression = _expression[:, index]
    if np.ptp(target_expression) == 0:
        return None
    with threadpool_limits(1):
        return core.infer_partial_network(
            regressor_type="GBM",
            regressor_kwargs=core.SGBM_KWARGS,
            tf_matrix=_tf_matrix,
            tf_matrix_gene_names=_tf_names,
            target_gene_name=target,
            target_gene_expression=target_expression,
            include_meta=False,
            early_stop_window_length=core.EARLY_STOP_WINDOW_LENGTH,
            seed=seed,
        )


networks = []
groups = adata.obs[obs_key].astype(str)
for group in sorted(groups.unique()):
    mask = (groups == group).to_numpy()
    n_cells = int(mask.sum())
    if n_cells < min_cells:
        print(f"Skipping group {group}: {n_cells} cells is below the minimum of {min_cells}.")
        continue

    _expression = dense(adata[mask, genes].X).astype(np.float64)
    _tf_matrix, _tf_names = core.to_tf_matrix(_expression, list(genes), list(regulators))

    with multiprocessing.get_context("fork").Pool(cpus) as pool:
        results = [result for result in pool.map(infer_target, target_indices) if result is not None]

    if not results:
        print(f"No links inferred for group {group}.")
        continue
    network = pd.concat(results, ignore_index=True)
    network = network.sort_values(["importance", "TF", "target"], ascending=[False, True, True])
    if max_links > 0:
        network = network.head(max_links)
    network.insert(0, "group", group)
    networks.append(network)
    print(f"Group {group}: {n_cells} cells, {len(network)} links kept.")

if networks:
    result = pd.concat(networks, ignore_index=True)
else:
    print(f"No group of {obs_key} reached {min_cells} cells; writing an empty network.")
    result = pd.DataFrame(
        {column: pd.Series(dtype="float64" if column == "importance" else "str") for column in COLUMNS}
    )

result = result[COLUMNS].astype({"group": str, "TF": str, "target": str, "importance": "float64"})
result.to_csv(f"{prefix}.tsv", sep="\\t", index=False)
result.reset_index(drop=True).reset_index(names="__index__").astype({"__index__": str}).to_parquet(
    Path(f"{prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "arboreto": version("arboreto"),
        "pandas": pd.__version__,
        "scikit-learn": version("scikit-learn"),
        "scanpy": version("scanpy"),
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
