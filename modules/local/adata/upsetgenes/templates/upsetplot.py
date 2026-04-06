#!/opt/conda/bin/python

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import base64
import json

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import matplotlib.pyplot as plt
import upsetplot
import matplotlib
import yaml

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)

prefix = "${prefix}"

def load_gene_names(path: str) -> list[str]:
    adata = ad.read_h5ad(path, backed="r")
    try:
        return adata.var_names.unique().tolist()
    finally:
        adata.file.close()

adata_genes = dict(zip(
    "${names.join(' ')}".split(),
    [load_gene_names(path) for path in "${h5ads}".split()]
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
