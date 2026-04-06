#!/opt/conda/bin/python

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import base64
import json

os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot
import matplotlib
import yaml

adata = ad.read_h5ad("${h5ad}")
threshold = int("${threshold}")
prefix = "${prefix}"

def load(path: str) -> pd.DataFrame:
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)

predictions = pd.concat([load(f) for f in "${predictions}".split()], axis=1)
mask = predictions.sum(axis=1) >= threshold

adata = adata[~mask, :]
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)

if not len(predictions.columns) > 1:
    exit(0)

# Plot

contents = {column: predictions[column][predictions[column]].index.tolist()
               for column in predictions.columns}

plot_data = upsetplot.from_contents(contents)

upsetplot.plot(plot_data,
               sort_by="cardinality",
               show_counts=True,
               subset_size="count",
               min_subset_size=10)
plot_path = f"{prefix}_predictions_mqc.png"
plt.savefig(plot_path)

# MultiQC

with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "${prefix}",
        "parent_id": "doublet_predictions",
        "parent_name": "Doublet predictions",
        "parent_description": "Upset plots of the various doublet prediction tools for each sample.",

        "section_name": "${meta.id}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)
