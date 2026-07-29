#!/usr/bin/env python3

import os
import platform
import warnings

import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

n_perms_raw = int("${n_perms}")
max_cells_raw = "${max_cells}"
subsample_strategy = "${subsample_strategy}"
subsample_seed = int("${subsample_seed}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"


def _optional_int(value):
    return None if value in ("", "null", "None") else int(value)


def _subsample_info(n_before, n_after, strategy, seed, max_cells):
    return {
        "subsampled": n_after < n_before,
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "strategy": strategy,
        "seed": seed,
        "max_cells": max_cells,
    }


def _subsample_groupby(strategy, obs_key, adata, batch_key="batch"):
    strategies = {
        "stratified_obs": [obs_key],
        "stratified_obs_batch": [obs_key, batch_key],
    }

    if strategy in ("none", ""):
        warnings.warn(
            "liana_max_cells is set but liana_subsample_strategy is 'none'; "
            "using stratified_obs instead of uniform sampling."
        )
        strategy = "stratified_obs"

    if strategy not in strategies:
        raise SystemExit(
            f"Unknown liana_subsample_strategy '{strategy}'; expected stratified_obs or stratified_obs_batch."
        )

    groupby_keys = strategies[strategy]
    if batch_key in groupby_keys and batch_key not in adata.obs:
        raise SystemExit(
            f"liana_subsample_strategy '{strategy}' requires obs column '{batch_key}'; "
            f"available: {list(adata.obs.columns)}"
        )
    if obs_key not in adata.obs:
        raise SystemExit(f"LIANA grouping column '{obs_key}' not in obs; available: {list(adata.obs.columns)}")

    return strategy, groupby_keys


def _stratified_subsample(adata, n_max, strategy, seed, obs_key):
    n_before = adata.n_obs
    info = _subsample_info(
        n_before=n_before,
        n_after=n_before,
        strategy="none",
        seed=seed,
        max_cells=n_max,
    )
    if n_max is None or n_max <= 0 or n_before <= n_max:
        return adata, info

    effective_strategy, groupby_keys = _subsample_groupby(strategy, obs_key, adata)
    rng = np.random.default_rng(seed)
    frac = n_max / n_before
    selected = np.concatenate(
        [
            rng.choice(idx, size=max(1, round(len(idx) * frac)), replace=False)
            for idx in adata.obs.groupby(groupby_keys, observed=True).indices.values()
        ]
    )

    if selected.size > n_max:
        selected = rng.choice(selected, size=n_max, replace=False)

    adata_sub = adata[np.sort(selected)].copy()
    return adata_sub, _subsample_info(
        n_before=n_before,
        n_after=adata_sub.n_obs,
        strategy=effective_strategy,
        seed=seed,
        max_cells=n_max,
    )


def _log_subsample_info(info):
    if info["subsampled"]:
        print(
            "LIANA subsampled "
            f"{info['n_cells_before']} -> {info['n_cells_after']} cells "
            f"(strategy={info['strategy']}, seed={info['seed']}, max_cells={info['max_cells']})."
        )
    else:
        print(f"LIANA using all {info['n_cells_before']} cells (no subsampling).")


max_cells = _optional_int(max_cells_raw)
n_perms = None if n_perms_raw <= 0 else n_perms_raw

if adata.obs[obs_key].nunique() > 1:
    adata, subsample_info = _stratified_subsample(
        adata,
        max_cells,
        subsample_strategy,
        subsample_seed,
        obs_key,
    )
    _log_subsample_info(subsample_info)

    sc.pp.log1p(adata)
    try:
        li.mt.rank_aggregate(
            adata,
            obs_key,
            use_raw=False,
            verbose=True,
            n_jobs=int("${task.cpus}"),
            n_perms=n_perms,
            seed=subsample_seed,
        )
        df: pd.DataFrame = adata.uns["liana_res"]

        df.to_pickle(f"{prefix}.pkl")
        adata.write_h5ad(f"{prefix}.h5ad")

    except ValueError as e:
        if "cannot set a frame with no defined index and a scalar" in str(e):
            print(f"Error: {e}")
        else:
            raise e
else:
    print(f"Skipping rank aggregation because the column {obs_key} has only one unique value.")

# Versions

versions = {
    "python": platform.python_version(),
    "scanpy": sc.__version__,
    "liana": li.__version__,
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
