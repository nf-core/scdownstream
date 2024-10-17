#!/opt/conda/bin/python

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot
import matplotlib
import platform
import base64
import json

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

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
    f.write(format_yaml_like(versions))

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
