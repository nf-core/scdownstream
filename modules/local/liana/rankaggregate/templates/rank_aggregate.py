#!/usr/bin/env python3

import base64
import json
import os
import platform
import warnings

import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import liana as li
import matplotlib.pyplot as plt
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


def write_mqc_plot(plot_path, plot_id, plot_label, section_name, description):
    with open(plot_path, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
    custom_json = {
        "id": plot_id,
        "parent_id": "${meta.integration}",
        "parent_name": "${meta.integration}",
        "parent_description": "Results of the ${meta.integration} integration.",
        "section_name": f"{section_name} ({plot_label})",
        "description": f"{description} {plot_label.capitalize()}.",
        "plot_type": "image",
        "data": image_html,
    }
    with open(f"{plot_id}_mqc.json", "w") as f_json:
        json.dump(custom_json, f_json)


def _sort_cell_labels(labels):
    def _key(label):
        try:
            return (0, float(label))
        except (TypeError, ValueError):
            return (1, str(label))

    return sorted((str(label) for label in labels), key=_key)


def _top_cell_types(df, n=6):
    counts = pd.concat([df["source"].astype(str), df["target"].astype(str)]).value_counts()
    return _sort_cell_labels(counts.head(n).index)


def _write_liana_plots(adata, df, obs_key):
    if df.empty:
        print("Skipping LIANA plots: liana_res is empty.")
        return

    section_name = f"LIANA (grouped by: {obs_key})"
    description = f"LIANA rank aggregate results, grouped by <code>{obs_key}</code>."
    has_specificity = "specificity_rank" in df.columns

    plot_df = df.copy()
    plot_df["source"] = plot_df["source"].astype(str)
    plot_df["target"] = plot_df["target"].astype(str)

    if has_specificity:
        filtered_df = plot_df.loc[plot_df["specificity_rank"] <= 0.05]
        if filtered_df.empty:
            filtered_df = plot_df
    else:
        filtered_df = plot_df

    # LIANA's faceted dotplot becomes unreadable with many groups; keep the most
    # interactive cell types (same pattern as the LIANA docs examples).
    n_groups = plot_df[["source", "target"]].stack().nunique()
    cell_types = _top_cell_types(filtered_df, n=6 if n_groups > 6 else n_groups)

    try:
        n_types = max(len(cell_types), 1)
        dotplot_kwargs = {
            "adata": adata,
            "liana_res": plot_df,
            "colour": "magnitude_rank",
            "inverse_colour": True,
            "top_n": 15,
            "orderby": "magnitude_rank",
            "orderby_ascending": True,
            "source_labels": cell_types,
            "target_labels": cell_types,
            "figure_size": (max(8.0, n_types * 1.8), max(6.0, 15 * 0.35 + 2)),
            "return_fig": True,
        }
        if has_specificity:
            dotplot_kwargs["size"] = "specificity_rank"
            dotplot_kwargs["inverse_size"] = True
            dotplot_kwargs["filter_fun"] = lambda x: x["specificity_rank"] <= 0.05
        else:
            dotplot_kwargs["size"] = "magnitude_rank"
            dotplot_kwargs["inverse_size"] = True

        dotplot_path = f"{prefix}_dotplot.png"
        fig = li.pl.dotplot(**dotplot_kwargs)
        fig.save(dotplot_path, dpi=150, verbose=False)
        write_mqc_plot(dotplot_path, f"{prefix}_dotplot", "dot plot", section_name, description)
    except Exception as e:
        print(f"Warning: LIANA dotplot failed: {e}")

    try:
        plot_adata = adata
        if not pd.api.types.is_string_dtype(adata.obs[obs_key]):
            plot_adata = adata.copy()
            plot_adata.obs[obs_key] = plot_adata.obs[obs_key].astype(str)

        circle_kwargs = {
            "adata": plot_adata,
            "liana_res": plot_df,
            "groupby": obs_key,
            "score_key": "magnitude_rank",
            "inverse_score": True,
            "pivot_mode": "counts",
            "figure_size": (10, 10),
            # Fixed cartesian offset drifts labels away from nodes on a circle.
            "node_label_offset": (0.0, 0.0),
            "node_label_size": 10,
        }
        if has_specificity:
            circle_kwargs["filter_fun"] = lambda x: x["specificity_rank"] <= 0.05

        ax = li.pl.circle_plot(**circle_kwargs)
        circle_path = f"{prefix}_circle.png"
        ax.figure.savefig(circle_path, dpi=150, bbox_inches="tight", pad_inches=0.3)
        plt.close(ax.figure)
        write_mqc_plot(circle_path, f"{prefix}_circle", "circle plot", section_name, description)
    except Exception as e:
        print(f"Warning: LIANA circle plot failed: {e}")

    try:
        tileplot_path = f"{prefix}_tileplot.png"
        plot = li.pl.tileplot(
            adata=adata,
            fill="means",
            label="props",
            label_fun=lambda x: f"{x:.2f}",
            top_n=10,
            orderby="magnitude_rank",
            orderby_ascending=True,
            uns_key="liana_res",
            figure_size=(8, 7),
        )
        plot.save(tileplot_path, dpi=150, verbose=False)
        write_mqc_plot(tileplot_path, f"{prefix}_tileplot", "tile plot", section_name, description)
    except Exception as e:
        print(f"Warning: LIANA tileplot failed: {e}")


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

    if (adata.X < 0).nnz == 0:
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

        df.to_parquet(f"{prefix}.parquet", index=True)
        adata.write_h5ad(f"{prefix}.h5ad")
        _write_liana_plots(adata, df, obs_key)

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
