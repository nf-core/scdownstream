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
context_key = "${context_key}"

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
            "LIANA by-sample subsampled "
            f"{info['n_cells_before']} -> {info['n_cells_after']} cells "
            f"(strategy={info['strategy']}, seed={info['seed']}, max_cells={info['max_cells']})."
        )
    else:
        print(f"LIANA by-sample using all {info['n_cells_before']} cells (no subsampling).")


def _build_contexts(adata, context_key):
    if context_key not in adata.obs.columns:
        raise SystemExit(f"Context column '{context_key}' not in obs; available: {list(adata.obs.columns)}")

    context_series = adata.obs[context_key].astype(str)
    if "condition" in adata.obs.columns:
        condition_series = adata.obs["condition"].astype(str)
        n_conditions = (
            pd.DataFrame({context_key: context_series, "condition": condition_series})
            .groupby(context_key, observed=True)["condition"]
            .nunique()
        )
        multi = n_conditions[n_conditions > 1]
        if not multi.empty:
            print(
                "Warning: omitting condition from context metadata because some contexts "
                f"map to multiple conditions: {multi.index.tolist()}"
            )
            contexts = pd.DataFrame({context_key: sorted(context_series.unique())})
        else:
            contexts = (
                pd.DataFrame({context_key: context_series, "condition": condition_series})
                .drop_duplicates()
                .sort_values(context_key)
                .reset_index(drop=True)
            )
    else:
        contexts = pd.DataFrame({context_key: sorted(context_series.unique())})

    return contexts


max_cells = _optional_int(max_cells_raw)
n_perms = None if n_perms_raw <= 0 else n_perms_raw

if context_key not in adata.obs.columns:
    raise SystemExit(f"Context column '{context_key}' not in obs; available: {list(adata.obs.columns)}")
if obs_key not in adata.obs.columns:
    raise SystemExit(f"LIANA grouping column '{obs_key}' not in obs; available: {list(adata.obs.columns)}")

n_contexts = adata.obs[context_key].nunique()
n_groups = adata.obs[obs_key].nunique()

if n_contexts < 2:
    print(
        f"Skipping LIANA by-sample because context column '{context_key}' "
        f"has {n_contexts} unique value(s); at least 2 are required."
    )
elif n_groups < 2:
    print(
        f"Skipping LIANA by-sample because grouping column '{obs_key}' "
        f"has {n_groups} unique value(s); at least 2 are required."
    )
else:
    adata, subsample_info = _stratified_subsample(
        adata,
        max_cells,
        subsample_strategy,
        subsample_seed,
        obs_key,
    )
    _log_subsample_info(subsample_info)

    sc.pp.log1p(adata)

    contexts = _build_contexts(adata, context_key)

    try:
        li.mt.rank_aggregate.by_sample(
            adata,
            groupby=obs_key,
            sample_key=context_key,
            use_raw=False,
            verbose=True,
            n_jobs=int("${task.cpus}"),
            n_perms=n_perms,
            seed=subsample_seed,
            inplace=True,
        )
        df: pd.DataFrame = adata.uns["liana_res"]
        if context_key not in df.columns:
            raise SystemExit(
                f"Expected context column '{context_key}' in by-sample results; columns: {list(df.columns)}"
            )

        df.to_csv(f"{prefix}.csv.gz", index=False, compression="gzip")
        contexts.to_csv(f"{prefix}_contexts.tsv", sep="\t", index=False)
    except ValueError as e:
        if "cannot set a frame with no defined index and a scalar" in str(e):
            print(f"Error: {e}")
        else:
            raise e

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "liana": li.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
