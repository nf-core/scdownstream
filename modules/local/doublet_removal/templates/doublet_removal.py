#!/opt/conda/bin/python

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot
import matplotlib
import platform

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

predictions = pd.concat([pd.read_pickle(f) for f in "${predictions}".split()], axis=1)
mask = predictions.sum(axis=1) >= threshold

adata = adata[mask, :]
adata.write_h5ad(f"{prefix}.h5ad")

print(predictions.head())

# Plot
plot_data = upsetplot.from_indicators(predictions)

upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, min_subset_size=10)
plot_path = f"{prefix}_predictions_mqc.png"
plt.savefig(plot_path)

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
