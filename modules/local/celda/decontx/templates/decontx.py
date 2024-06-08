#!/usr/bin/env python3

import anndata as ad
import anndata2ri
import rpy2
import rpy2.robjects as ro
import platform
celda = ro.packages.importr('celda')

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
sc_experiment = anndata2ri.py2rpy(adata)

corrected = celda.decontX(sc_experiment, batch=adata.obs['${batch_col}'].tolist())
counts = celda.decontXcounts(corrected)

adata.layers['ambient'] = anndata2ri.rpy2py(counts).T
adata.X = adata.layers['ambient']
adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": rpy2.__version__,
        "celda": celda.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
