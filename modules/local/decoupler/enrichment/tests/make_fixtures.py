"""Generate small edgePython-style DE tables (no stat column) for the enrichment tests."""

from pathlib import Path

import numpy as np
import pandas as pd

outdir = Path(__file__).parent
network = pd.read_csv(outdir.parents[1] / "network" / "tests" / "network.tsv", sep="\t")
genes = sorted(set(network["target"]) - {"NOTAGENE"}) + [f"GENE{i}" for i in range(40)]
rng = np.random.default_rng(42)

for celltype, shift in [("B_cell", 1.5), ("Monocyte", -1.0)]:
    log2fc = rng.normal(0, 1, len(genes))
    log2fc[[genes.index(g) for g in ("CD79A", "IGHM", "CD22", "FCRL1", "CR2")]] += shift
    pvalue = np.clip(np.exp(-np.abs(log2fc) * 3), 1e-12, 1)
    table = pd.DataFrame(
        {
            "gene": genes,
            "log2fc": log2fc.round(4),
            "pvalue": pvalue.round(8),
            "padj": np.clip(pvalue * 2, 0, 1).round(8),
            "group": "Disease",
            "contrast": "Disease vs Healthy",
            "stratum": f"celltype={celltype}",
        }
    )
    table.to_parquet(outdir / f"test_{celltype}_Disease_results.parquet", index=False)
