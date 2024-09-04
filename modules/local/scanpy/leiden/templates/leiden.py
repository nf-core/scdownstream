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
use_gpu = "${task.ext.use_gpu}" == "true"
prefix = "${prefix}"

kwargs = {
    "resolution": float("${resolution}"),
    "key_added": prefix
}

if use_gpu:
    import rapids_singlecell as rsc
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

    rsc.get.anndata_to_GPU(adata)
    rsc.tl.leiden(adata, **kwargs)
    rsc.get.anndata_to_CPU(adata)
else:
    sc.tl.leiden(adata, **kwargs)

adata.obs[[prefix]].to_pickle(f"{prefix}.pkl")
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
