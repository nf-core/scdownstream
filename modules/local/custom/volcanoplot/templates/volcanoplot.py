#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import base64
import json
import platform
import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import yaml

# Glob the staged parquet so Nextflow does not escape colons in filenames
# such as condition:leiden:0:wilcoxon_characteristic_genes_results.parquet.
staged_parquets = sorted(Path(".").glob("*.parquet"))
if not staged_parquets:
    raise FileNotFoundError("No staged DE parquet file found")
parquet_path = staged_parquets[0].resolve()
prefix = "${prefix}"
interesting_path = "${interesting_genes}"
integration = "${meta.integration}"
de_method = "${meta.de_method}"

REQUIRED_COLUMNS = ["gene", "log2fc", "pvalue", "padj", "group", "contrast", "stratum"]
SCANPY_METHODS = {"wilcoxon", "t-test", "t-test_overestim_var", "logreg"}
METHOD_LABELS = {
    "pydeseq2": "PyDESeq2",
    "edgepython": "edgePython",
    "edgepython_sc": "edgepython_sc",
}


def load_interesting_genes(path_str):
    if not path_str or path_str in ("[]", "null", "None"):
        return set()
    genes = set()
    first_data = True
    with open(path_str) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            field = line.split(",")[0].strip().strip('"').strip("'")
            if not field:
                continue
            if first_data and field.lower() in ("gene", "symbol", "names"):
                first_data = False
                continue
            first_data = False
            genes.add(field.lower())
    return genes


def method_label_for(method):
    if not method or method in ("null", "None"):
        return "DE"
    if method in SCANPY_METHODS:
        return f"Scanpy {method}"
    return METHOD_LABELS.get(method, method)


def method_slug_for(method):
    if method in SCANPY_METHODS:
        slug = f"scanpy_{method}"
    elif method and method not in ("null", "None"):
        slug = method
    else:
        slug = "de"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", slug)


def mqc_parent(method_slug, method_label):
    parent_id = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{integration}_{method_slug}")
    return {
        "parent_id": parent_id,
        "parent_name": f"{integration}: {method_label}",
        "parent_description": f"Differential expression volcano plots from {method_label} ({integration} integration).",
    }


def unique_nonempty(series):
    values = []
    seen = set()
    for value in series.tolist():
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none", "null"}:
            continue
        if text not in seen:
            seen.add(text)
            values.append(text)
    return values


def single_stratum(frame):
    strata = unique_nonempty(frame["stratum"])
    return strata[0] if len(strata) == 1 else ""


def split_vs(contrast):
    marker = " vs "
    if marker in contrast:
        left, right = contrast.split(marker, 1)
        return left, right
    return contrast, None


def output_stem(path):
    stem = path.stem
    suffix = "_results"
    if stem.endswith(suffix):
        return f"{stem[: -len(suffix)]}_volcano"
    if stem:
        return f"{stem}_volcano"
    return f"{prefix}_volcano"


def prepare_volcano_df(df, interesting):
    plot_df = df.copy()
    if "gene" not in plot_df.columns or "log2fc" not in plot_df.columns:
        return None, None, None

    padj_vals = pd.to_numeric(plot_df["padj"], errors="coerce") if "padj" in plot_df.columns else None
    pvalue_vals = pd.to_numeric(plot_df["pvalue"], errors="coerce") if "pvalue" in plot_df.columns else None
    use_padj = padj_vals is not None and padj_vals.notna().any()
    if use_padj:
        p_use = "padj"
        plot_df["padj"] = padj_vals
    elif pvalue_vals is not None:
        p_use = "pvalue"
        plot_df["pvalue"] = pvalue_vals
    else:
        return None, None, None

    plot_df["gene"] = plot_df["gene"].astype(str)
    plot_df["log2fc"] = pd.to_numeric(plot_df["log2fc"], errors="coerce")
    plot_df[p_use] = pd.to_numeric(plot_df[p_use], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        plot_df["neglog10"] = -np.log10(plot_df[p_use])
    plot_df["neglog10"] = plot_df["neglog10"].replace([np.inf, -np.inf], np.nan)
    plot_df = plot_df.dropna(subset=["log2fc", "neglog10"])
    if len(plot_df) < 2:
        return None, None, None

    if use_padj:
        padj_plot = pd.to_numeric(plot_df["padj"], errors="coerce")
        significant = (padj_plot < 0.05) & (plot_df["log2fc"].abs() >= 1)
    else:
        significant = (plot_df[p_use] < 0.05) & (plot_df["log2fc"].abs() >= 1)
    plot_df["significant"] = significant.fillna(False)
    plot_df["interesting"] = plot_df["gene"].str.lower().isin(interesting)
    return plot_df, p_use, use_padj


def _genes_to_label(plot_df, p_use, max_labels=10):
    """Label significant genes only. Prefer interesting when provided; never label grey points."""
    significant = plot_df[plot_df["significant"]].sort_values(p_use, ascending=True)
    if significant.empty:
        return significant
    if plot_df["interesting"].any():
        return significant[significant["interesting"]].head(max_labels)
    return significant.head(max_labels)


def _annotate_genes(ax, to_label, fontsize=6):
    if to_label.empty:
        return
    texts = [
        ax.text(row["log2fc"], row["neglog10"], row["gene"], fontsize=fontsize, alpha=0.9)
        for _, row in to_label.iterrows()
    ]
    try:
        from adjustText import adjust_text

        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.35),
            expand=(1.2, 1.4),
            force_text=(0.5, 0.8),
            ensure_inside_axes=True,
        )
    except Exception:
        pass


def draw_volcano_ax(ax, plot_df, p_use, use_padj, title, marker_s=(10, 28), xlabel="log2FC", ylabel=None, label_size=6):
    colours = plot_df["significant"].map({True: "#d62728", False: "#7f7f7f"})
    other = plot_df[~plot_df["interesting"]]
    marked = plot_df[plot_df["interesting"]]
    if not other.empty:
        ax.scatter(
            other["log2fc"],
            other["neglog10"],
            c=colours.loc[other.index],
            marker="o",
            alpha=0.5,
            s=marker_s[0],
            linewidths=0,
        )
    if not marked.empty:
        ax.scatter(
            marked["log2fc"],
            marked["neglog10"],
            c=colours.loc[marked.index],
            marker="^",
            alpha=0.9,
            s=marker_s[1],
            linewidths=0,
        )
    ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.axvline(-1, color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.axvline(1, color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax.set_xlabel(xlabel, fontsize=8)
    if ylabel is None:
        ylabel = "-log10(padj)" if use_padj else "-log10(p)"
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    _annotate_genes(ax, _genes_to_label(plot_df, p_use), fontsize=label_size)


def write_mqc_image(out_png, out_mqc_id, section_name, description, parent):
    with open(out_png, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
    custom_json = {
        "id": out_mqc_id,
        "parent_id": parent["parent_id"],
        "parent_name": parent["parent_name"],
        "parent_description": parent["parent_description"],
        "section_name": section_name,
        "description": description,
        "plot_type": "image",
        "data": image_html,
    }
    with open(f"{out_mqc_id}_mqc.json", "w") as f_json:
        json.dump(custom_json, f_json)


def write_volcano(plot_df, p_use, use_padj, out_png, out_mqc_id, section_name, description, parent):
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ylabel = "-log10(adjusted p-value)" if use_padj else "-log10(p-value)"
    colours = plot_df["significant"].map({True: "#d62728", False: "#7f7f7f"})
    other = plot_df[~plot_df["interesting"]]
    marked = plot_df[plot_df["interesting"]]
    if not other.empty:
        ax.scatter(
            other["log2fc"], other["neglog10"], c=colours.loc[other.index], marker="o", alpha=0.5, s=16, linewidths=0
        )
    if not marked.empty:
        ax.scatter(
            marked["log2fc"],
            marked["neglog10"],
            c=colours.loc[marked.index],
            marker="^",
            alpha=0.9,
            s=36,
            linewidths=0,
        )
    ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.axvline(-1, color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel(ylabel)
    ax.set_title(section_name)
    _annotate_genes(ax, _genes_to_label(plot_df, p_use), fontsize=7)
    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    write_mqc_image(out_png, out_mqc_id, section_name, description, parent)


def write_volcano_grid(panels, out_png, out_mqc_id, section_name, description, parent):
    n = len(panels)
    ncols = min(4, max(1, int(np.ceil(np.sqrt(n)))))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.6 * ncols, 3.0 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    for idx, (title, plot_df, p_use, use_padj) in enumerate(panels):
        row, col = divmod(idx, ncols)
        draw_volcano_ax(axes[row][col], plot_df, p_use, use_padj, title)
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    fig.suptitle(section_name, fontsize=11)
    plt.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    write_mqc_image(out_png, out_mqc_id, section_name, description, parent)


def section_with_stratum(base, stratum):
    return f"{base} (within {stratum})" if stratum else base


def description_contrast(method_label, left, right, stratum):
    if right:
        text = f"{method_label} volcano for contrast <code>{left}</code> versus <code>{right}</code>"
    else:
        text = f"{method_label} volcano for contrast <code>{left}</code>"
    if stratum:
        text += f" within <code>{stratum}</code>"
    return f"{text}."


def description_grid(method_label, stratum):
    text = f"{method_label} volcano panels for each group versus the rest of the cells"
    if stratum:
        text += f" within <code>{stratum}</code>"
    return f"{text}."


out_stem = output_stem(parquet_path)
out_png = f"{out_stem}.png"
out_mqc_id = out_stem

df = pd.read_parquet(parquet_path)
missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
if missing:
    raise ValueError(f"Parquet is missing required columns: {', '.join(missing)}")

interesting = load_interesting_genes(interesting_path)
method_label = method_label_for(de_method)
volcano_parent = mqc_parent(method_slug_for(de_method), method_label)

padj_all_null = pd.to_numeric(df["padj"], errors="coerce").notna().sum() == 0
pvalue_all_null = pd.to_numeric(df["pvalue"], errors="coerce").notna().sum() == 0
if padj_all_null and pvalue_all_null:
    print(f"Warning: missing p-values for {out_png}; skipping volcano.")
else:
    panels = []
    for group_name, gdf in df.groupby("group", sort=True):
        plot_df, p_use, use_padj = prepare_volcano_df(gdf, interesting)
        if plot_df is None:
            continue
        panels.append((str(group_name), plot_df, p_use, use_padj))

    if not panels:
        print(f"Warning: no plottable volcano panels for {out_png}; skipping.")
    else:
        stratum = single_stratum(df)
        if len(panels) == 1:
            _group_name, plot_df, p_use, use_padj = panels[0]
            contrasts = unique_nonempty(df["contrast"])
            contrast = contrasts[0] if contrasts else str(_group_name)
            section_name = section_with_stratum(f"{method_label} volcano: {contrast}", stratum)
            left, right = split_vs(contrast)
            description = description_contrast(method_label, left, right, stratum)
            write_volcano(
                plot_df,
                p_use,
                use_padj,
                out_png,
                out_mqc_id,
                section_name,
                description,
                volcano_parent,
            )
        elif len(panels) == 2:
            # Two groups: A vs rest and B vs rest are mirrors; keep a single A vs B panel.
            group_a, plot_df, p_use, use_padj = panels[0]
            group_b = panels[1][0]
            contrast = f"{group_a} vs {group_b}"
            section_name = section_with_stratum(f"{method_label} volcano: {contrast}", stratum)
            description = description_contrast(method_label, group_a, group_b, stratum)
            write_volcano_grid(
                [(contrast, plot_df, p_use, use_padj)],
                out_png,
                out_mqc_id,
                section_name,
                description,
                volcano_parent,
            )
        else:
            grid_panels = [
                (f"{group_name} vs rest", plot_df, p_use, use_padj) for group_name, plot_df, p_use, use_padj in panels
            ]
            contrasts = unique_nonempty(df["contrast"])
            contrast_label = contrasts[0] if len(contrasts) == 1 else "vs rest"
            section_name = section_with_stratum(f"{method_label} volcanoes: {contrast_label}", stratum)
            description = description_grid(method_label, stratum)
            write_volcano_grid(
                grid_panels,
                out_png,
                out_mqc_id,
                section_name,
                description,
                volcano_parent,
            )

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "pyarrow": pa.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
