#!/usr/bin/env python

import os
import platform
import base64
import json

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import matplotlib.pyplot as plt
import upsetplot
import matplotlib

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

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))

prefix = "${prefix}"

adata_genes = dict(zip(
    "${names.join(' ')}".split(),
    [sc.read_h5ad(path, backed='r').var.index.unique().to_list() for path in "${h5ads}".split()]
))

if len(adata_genes) < 2:
    exit(0)

plot_data = upsetplot.from_contents(adata_genes)

upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, min_subset_size=10)
plot_path = f"{prefix}_genes.png"
plt.savefig(plot_path)

# MultiQC

with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "${prefix}",
        "section_name": "Genes upset: ${prefix}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)
