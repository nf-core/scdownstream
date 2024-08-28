#!/usr/bin/env python3

import scanpy as sc
import platform

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

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
prefix = "${prefix}"
n_hvgs = int("${n_hvgs}")

adata.layers["log1p"] = sc.pp.log1p(adata.X)

if adata.n_vars > n_hvgs:
    sc.pp.highly_variable_genes(adata,
                                n_top_genes=int("${n_hvgs}"),
                                batch_key="batch",
                                layer="log1p")
    adata = adata[:, adata.var["highly_variable"]]

del adata.layers["log1p"]

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
