#!/usr/bin/env python3

import os
from importlib.metadata import version

os.environ["TMPDIR"] = "."

import anndata as ad
import anndata2ri
import pandas as pd
import rpy2.robjects as ro
import yaml

seurat = ro.packages.importr("Seurat")
# Import SingleCellExperiment to check for class
sce_pkg = ro.packages.importr("SingleCellExperiment")

# Read the RDS file first
rds_obj = ro.r('readRDS("${rds}")')

# Check if it's already a SingleCellExperiment
is_sce = "SingleCellExperiment" in rds_obj.extends()

# Convert only if not already a SingleCellExperiment
if not is_sce:
    sce = seurat.as_SingleCellExperiment(rds_obj)
else:
    sce = rds_obj

adata = anndata2ri.rpy2py(sce)

# Convert indices to string
adata.obs.index = adata.obs.index.astype(str)
adata.var.index = adata.var.index.astype(str)

adata.write_h5ad("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "anndata": ad.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": version("rpy2"),
        "pandas": pd.__version__,
        "seurat": seurat.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
