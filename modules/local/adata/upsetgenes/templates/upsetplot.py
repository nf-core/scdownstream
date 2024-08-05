#!/opt/conda/bin/python

import scanpy as sc
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

adata = sc.read_h5ad("${h5ad}")
split_col = "${split_col}"
prefix = "${prefix}"
sample_genes = {}

# Split into multiple adatas, based on sample
samples = adata.obs[split_col].unique()

if len(samples) < 2:
    exit(0)

for sample in samples:
    adata_sample = adata[adata.obs[split_col] == sample].copy()
    # Keep only genes with at least 1 count in at least 1 cell
    sc.pp.filter_genes(adata_sample, min_cells=1)
    sample_genes[sample] = set(adata_sample.var_names)

plot_data = upsetplot.from_contents(sample_genes)

upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, min_subset_size=10)
plot_path = f"{prefix}_{split_col}_genes.png"
plt.savefig(plot_path)

# MultiQC

with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "upset_${prefix}",
        "section_name": "Genes upset: ${prefix}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)
