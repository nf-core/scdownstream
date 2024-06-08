#!/usr/bin/env python3

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
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

adata = sc.read_h5ad("${h5ad}")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

print(adata)

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
plt.savefig("${prefix}_total_counts_vs_n_genes_by_counts_mqc.png")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
