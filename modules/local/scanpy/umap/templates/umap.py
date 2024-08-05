#!/usr/bin/env python3

import base64
import json
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
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

sc.tl.umap(adata)

adata.write_h5ad(f"{prefix}.h5ad")
df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Save UMAP as png, color by "batch"
plot_path = f"{prefix}_batch.png"
sc.pl.umap(adata, color="batch")
plt.savefig(plot_path)

# MultiQC

with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "${prefix}",
        "parent_id": "integration_${meta.integration}",
        "parent_name": "Integration: ${meta.integration}",
        "parent_description": "Plots illustrating the results of the ${meta.integration} integration.",
        "section_name": "UMAP (batch)",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
