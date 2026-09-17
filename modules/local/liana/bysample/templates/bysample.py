#!/usr/bin/env python3

import os
import platform
import warnings

import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import liana as li
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

# Match LIANA's default: cell identities with fewer cells are dropped per context.
LIANA_MIN_CELLS = 5


def _optional_int(value):
    return None if value in ("", "null", "None") else int(value)


def _eligible_contexts(obs, context_key, obs_key, min_cells):
    """Contexts with at least two cell groups that each have >= min_cells cells."""
    counts = (
        pd.DataFrame(
            {
                "context": obs[context_key].astype(str),
                "group": obs[obs_key].astype(str),
            }
        )
        .groupby(["context", "group"], observed=True)
        .size()
    )
    n_qualifying_groups = counts[counts >= min_cells].groupby("context").size()
    return set(n_qualifying_groups[n_qualifying_groups >= 2].index.tolist())


def _drop_ineligible_contexts(obs, context_key, obs_key, min_cells, stage):
    context_str = obs[context_key].astype(str)
    eligible = _eligible_contexts(obs, context_key, obs_key, min_cells)
    dropped = sorted(set(context_str.unique()) - eligible)
    if dropped:
        print(
            f"Dropping {len(dropped)} context(s) with too little data for by-sample LIANA "
            f"({stage}; need >=2 groups with >= {min_cells} cells each): {dropped}"
        )
    return obs.loc[context_str.isin(eligible)]


def _subsample_groupby_keys(strategy, obs_columns, obs_key, context_key, batch_key="batch"):
    # Always stratify by context so rare patients are not wiped out globally.
    strategies = {
        "stratified_obs": [context_key, obs_key],
        "stratified_obs_batch": [context_key, obs_key, batch_key],
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
    missing = [key for key in groupby_keys if key not in obs_columns]
    if missing:
        raise SystemExit(
            f"liana_subsample_strategy '{strategy}' requires obs column(s) {missing}; available: {list(obs_columns)}"
        )
    return strategy, groupby_keys


def _stratified_subsample_obs(obs, n_max, strategy, seed, obs_key, context_key):
    n_before = len(obs)
    if n_max is None or n_max <= 0 or n_before <= n_max:
        return obs, "none", n_before

    effective_strategy, groupby_keys = _subsample_groupby_keys(strategy, obs.columns, obs_key, context_key)
    sampled = obs.groupby(groupby_keys, observed=True, group_keys=False).sample(
        frac=n_max / n_before,
        random_state=seed,
    )
    if len(sampled) > n_max:
        sampled = sampled.sample(n=n_max, random_state=seed)
    return sampled, effective_strategy, n_before


def _build_contexts(obs, context_key):
    context_series = obs[context_key].astype(str)
    if "condition" in obs.columns:
        condition_series = obs["condition"].astype(str)
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
            return pd.DataFrame({context_key: sorted(context_series.unique())})
        return (
            pd.DataFrame({context_key: context_series, "condition": condition_series})
            .drop_duplicates()
            .sort_values(context_key)
            .reset_index(drop=True)
        )
    return pd.DataFrame({context_key: sorted(context_series.unique())})


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
    # Decide the final cell set using obs only, then materialise AnnData once.
    # Deduplicate so context_key == "condition" (or "batch") does not yield a
    # DataFrame on column select, which breaks pd.DataFrame({...}) construction.
    obs_cols = list(
        dict.fromkeys(
            [context_key, obs_key]
            + (["condition"] if "condition" in adata.obs.columns else [])
            + (["batch"] if "batch" in adata.obs.columns else [])
        )
    )
    obs = adata.obs.loc[:, obs_cols]

    obs = _drop_ineligible_contexts(obs, context_key, obs_key, LIANA_MIN_CELLS, stage="before subsampling")
    if obs[context_key].nunique() < 2:
        print(
            "Skipping LIANA by-sample because fewer than 2 contexts remain after dropping "
            f"contexts with too little data ({obs[context_key].nunique()} left)."
        )
    else:
        obs, strategy_used, n_before = _stratified_subsample_obs(
            obs,
            max_cells,
            subsample_strategy,
            subsample_seed,
            obs_key,
            context_key,
        )
        if len(obs) < n_before:
            print(
                "LIANA by-sample subsampled "
                f"{n_before} -> {len(obs)} cells "
                f"(strategy={strategy_used}, seed={subsample_seed}, max_cells={max_cells})."
            )
        else:
            print(f"LIANA by-sample using all {n_before} cells (no subsampling).")

        obs = _drop_ineligible_contexts(obs, context_key, obs_key, LIANA_MIN_CELLS, stage="after subsampling")
        if obs[context_key].nunique() < 2:
            print(
                "Skipping LIANA by-sample because fewer than 2 contexts remain after "
                f"subsampling and dropping sparse contexts ({obs[context_key].nunique()} left)."
            )
        else:
            contexts = _build_contexts(obs, context_key)
            keep_indexer = adata.obs_names.get_indexer(obs.index)
            del obs
            # In-place subset avoids holding the full object and a copy at once.
            adata._inplace_subset_obs(keep_indexer)

            sc.pp.log1p(adata)

            li.mt.rank_aggregate.by_sample(
                adata,
                groupby=obs_key,
                sample_key=context_key,
                use_raw=False,
                verbose=True,
                n_jobs=int("${task.cpus}"),
                n_perms=n_perms,
                seed=subsample_seed,
                min_cells=LIANA_MIN_CELLS,
                inplace=True,
            )
            df: pd.DataFrame = adata.uns["liana_res"]
            if context_key not in df.columns:
                raise SystemExit(
                    f"Expected context column '{context_key}' in by-sample results; columns: {list(df.columns)}"
                )

            df.to_csv(f"{prefix}.csv.gz", index=False, compression="gzip")
            contexts.to_csv(f"{prefix}_contexts.tsv", sep="\\t", index=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "liana": li.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
