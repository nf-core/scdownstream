#!/usr/bin/env python3

import base64
import json
import pickle
import platform
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

method = "${method}"
tmin = int("${tmin}")
prefix = "${prefix}"
de_method = "${meta.de_method ?: ''}"
TOP_SOURCES = 20

if method not in dc.mt.show()["name"].tolist():
    raise ValueError(f"Unknown decoupler method '{method}'")
run_method = getattr(dc.mt, method)


def gene_stat(table: pd.DataFrame) -> pd.Series:
    """Use the DE test statistic when present, otherwise the signed -log10 p-value."""
    if "stat" in table.columns and table["stat"].notna().any():
        stat = pd.to_numeric(table["stat"], errors="coerce")
    else:
        pvalue = pd.to_numeric(table["pvalue"], errors="coerce").clip(lower=np.finfo(float).tiny)
        stat = np.sign(pd.to_numeric(table["log2fc"], errors="coerce")) * -np.log10(pvalue)
    stat.index = table["gene"].astype(str)
    stat = stat.replace([np.inf, -np.inf], np.nan).dropna()
    return stat[~stat.index.duplicated()]


def first_value(table: pd.DataFrame, column: str) -> str:
    return str(table[column].iloc[0]) if column in table.columns and len(table) else ""


contrasts = []
for de_path in sorted(Path("de").glob("*.parquet")):
    table = pd.read_parquet(de_path)
    stat = gene_stat(table)
    if stat.empty:
        print(f"Warning: skipped {de_path.name} without usable statistics")
        continue
    info = {col: first_value(table, col) for col in ("stratum", "contrast", "group")}
    contrasts.append((info, stat.to_frame(de_path.stem).T))


def enrich(net: pd.DataFrame, resource: str) -> pd.DataFrame:
    rows = []
    for info, data in contrasts:
        try:
            scores, pvalues = run_method(data, net, tmin=tmin)
        except AssertionError as exc:
            print(f"Warning: skipped {resource} for {data.index[0]}: {exc}")
            continue
        padj = pvalues.iloc[0].to_numpy() if pvalues is not None else np.full(scores.shape[1], np.nan)
        if method == "mlm" and pvalues is not None:
            # decoupler adjusts p-values for all methods except mlm
            padj = false_discovery_control(padj)
        rows.append(pd.DataFrame({"source": scores.columns, "score": scores.iloc[0].to_numpy(), "padj": padj, **info}))
    if not rows:
        return pd.DataFrame()
    long = pd.concat(rows, ignore_index=True)
    long.insert(0, "resource", resource)
    long["de_method"] = de_method
    return long


def plot_dots(long: pd.DataFrame, resource: str) -> str:
    long = long.assign(label=long["stratum"].str.split("=", n=1).str[-1] + ": " + long["contrast"])
    top = (
        long.groupby("source")["score"].apply(lambda s: s.abs().max()).sort_values(ascending=False).index[:TOP_SOURCES]
    )
    sub = long[long["source"].isin(top)]
    sources, labels = sorted(top), sorted(sub["label"].unique())
    x = sub["label"].map({label: i for i, label in enumerate(labels)})
    y = sub["source"].map({source: i for i, source in enumerate(sources)})
    sizes = 20 + 30 * -np.log10(sub["padj"].fillna(1).clip(lower=1e-10))
    limit = float(sub["score"].abs().max()) or 1.0
    fig, ax = plt.subplots(
        figsize=(max(5, 0.6 * len(labels) + 4), max(4, 0.35 * len(sources) + 3)), constrained_layout=True
    )
    points = ax.scatter(x, y, c=sub["score"], s=sizes, cmap="RdBu_r", vmin=-limit, vmax=limit, edgecolors="grey")
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(sources)), sources)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(-0.5, len(sources) - 0.5)
    ax.set_title(f"{resource} enrichment ({method}); size: -log10 padj")
    fig.colorbar(points, ax=ax, label="Score")
    plot_path = f"{prefix}_{resource}_enrichment.png"
    fig.savefig(plot_path)
    plt.close(fig)
    return plot_path


def write_mqc(plot_path: str, resource: str) -> None:
    with open(plot_path, "rb") as f_plot, open(f"{prefix}_{resource}_enrichment_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        custom_json = {
            "id": f"{prefix}_{resource}_enrichment",
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Pathway and transcription factor enrichment of differential expression results for ${meta.integration}.",
            "section_name": " ".join(filter(None, ["${meta.id}", de_method, resource, "enrichment"])),
            "plot_type": "image",
            "data": f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>',
        }
        json.dump(custom_json, f_json)


results = {}
for network_path in sorted(Path("networks").glob("*.tsv")):
    resource = network_path.stem
    long = enrich(pd.read_csv(network_path, sep="\\t"), resource)
    if long.empty:
        continue
    results[resource] = long
    write_mqc(plot_dots(long, resource), resource)

if results:
    combined = pd.concat(results.values(), ignore_index=True)
    # Absolute path so pyarrow does not treat colons in Nextflow prefixes as URI schemes
    combined.to_parquet(str(Path(f"{prefix}_enrichment.parquet").resolve()), index=False, engine="pyarrow")

with open(f"{prefix}_decoupler_enrichment.pkl", "wb") as f:
    pickle.dump(results, f)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "decoupler": dc.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
