#!/usr/bin/env python3

import base64
import json
import pickle
import platform
from pathlib import Path

import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = ad.read_h5ad("${h5ad}")
method = "${method}"
tmin = int("${tmin}")
prefix = "${prefix}"
meta_id = "${meta.id}"
OBS_COLS = ["donor", "celltype", "condition", "n_cells"]
TOP_SOURCES = 25

if method not in dc.mt.show()["name"].tolist():
    raise ValueError(f"Unknown decoupler method '{method}'")
run_method = getattr(dc.mt, method)

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
data = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
obs = adata.obs[[col for col in OBS_COLS if col in adata.obs.columns]].copy()


def run_network(net: pd.DataFrame) -> pd.DataFrame:
    scores, pvalues = run_method(data, net, tmin=tmin)
    if pvalues is None:
        pvalues = pd.DataFrame(np.nan, index=scores.index, columns=scores.columns)
    elif method == "mlm":
        # decoupler adjusts p-values for all methods except mlm
        pvalues = pd.DataFrame(
            false_discovery_control(pvalues.values, axis=1), index=pvalues.index, columns=pvalues.columns
        )
    long = scores.rename_axis("sample_id").melt(ignore_index=False, var_name="source", value_name="score")
    long["padj"] = pvalues.melt(value_name="padj")["padj"].to_numpy()
    long = long.reset_index().merge(obs, left_on="sample_id", right_index=True, how="left")
    return long[["sample_id", *obs.columns, "source", "score", "padj"]]


def plot_heatmap(long: pd.DataFrame, resource: str) -> str:
    group_cols = [col for col in ("celltype", "condition") if col in long.columns] or ["sample_id"]
    means = long.pivot_table(index=group_cols, columns="source", values="score", aggfunc="mean", observed=True)
    top = means.var(axis=0).sort_values(ascending=False).index[:TOP_SOURCES]
    means = means[sorted(top)]
    limit = float(np.nanmax(np.abs(means.values))) or 1.0
    fig, ax = plt.subplots(
        figsize=(max(6, 0.35 * means.shape[1] + 3), max(4, 0.3 * means.shape[0] + 2)), constrained_layout=True
    )
    image = ax.imshow(means.values, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(means.shape[1]), means.columns, rotation=90)
    ax.set_yticks(
        range(means.shape[0]), [" | ".join(map(str, idx if isinstance(idx, tuple) else (idx,))) for idx in means.index]
    )
    ax.set_title(f"{resource} activity ({method})")
    fig.colorbar(image, ax=ax, label="Mean score")
    plot_path = f"{prefix}_{resource}_activity.png"
    fig.savefig(plot_path)
    plt.close(fig)
    return plot_path


def write_mqc(plot_path: str, resource: str) -> None:
    with open(plot_path, "rb") as f_plot, open(f"{prefix}_{resource}_activity_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        custom_json = {
            "id": f"{prefix}_{resource}_activity",
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Pathway and transcription factor activity for ${meta.integration}.",
            "section_name": f"{meta_id} {resource} activity",
            "plot_type": "image",
            "data": f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>',
        }
        json.dump(custom_json, f_json)


results = {}
for network_path in sorted(Path("networks").glob("*.tsv")):
    resource = network_path.stem
    net = pd.read_csv(network_path, sep="\\t")
    try:
        long = run_network(net)
    except AssertionError as exc:
        print(f"Warning: skipped network {resource}: {exc}")
        continue
    long.insert(0, "resource", resource)
    results[resource] = long
    # Absolute path so pyarrow does not treat colons in Nextflow prefixes as URI schemes
    long.to_parquet(str(Path(f"{prefix}_{resource}_scores.parquet").resolve()), index=False, engine="pyarrow")
    write_mqc(plot_heatmap(long, resource), resource)

with open(f"{prefix}_decoupler_activity.pkl", "wb") as f:
    pickle.dump(results, f)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "decoupler": dc.__version__,
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
