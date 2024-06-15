#!/usr/bin/env python3

import scanpy
import anndata2ri
import rpy2
import rpy2.robjects as ro
import platform
soupx = ro.packages.importr('soupx')

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

adata = scanpy.read_h5ad("${h5ad}")

adata_pp = adata.copy()
scanpy.pp.normalize_per_cell(adata_pp)
scanpy.pp.log1p(adata_pp)

scanpy.pp.pca(adata_pp)
scanpy.pp.neighbors(adata_pp)
scanpy.tl.leiden(adata_pp, key_added="soupx_groups")

soupx_groups = adata_pp.obs["soupx_groups"]
del adata_pp

cells = adata.obs_names
genes = adata.var_names
data = adata.X.T

adata_raw = scanpy.read_h5ad("${raw}")
data_raw = adata_raw.X.T
del adata_raw

r_data = anndata2ri.py2rpy(data)
r_data_raw = anndata2ri.py2rpy(data_raw)
ro.r("function(data, genes, cells) { rownames(data) <- genes; colnames(data) <- cells }")(data, genes, cells)

sc = soupx.SoupChannel(r_data_raw, r_data, calcSoupProfile = False)

soupProf = ro.r("function(data) { data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data)) }")(r_data)
sc = soupx.setSoupProfile(sc, soupProf)
sc = soupx.setClusters(sc, soupx_groups)
sc = soupx.autoEstCont(sc, doPlot = False)
out = soupx.adjustCounts(sc, roundToInt = True)

adata.layers["ambient"] = out.T

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": scanpy.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": rpy2.__version__,
        "soupx": soupx.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
