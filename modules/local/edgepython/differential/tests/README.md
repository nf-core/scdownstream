# Dispersion-design regression

`design.nf.test` executes the production pseudobulk and differential-expression processes. The paired case uses three donors, each represented in both conditions. The unpaired case adds `unpaired_donor` through `ADATA_EXTEND`, then reaggregates the existing pseudobulks under condition-specific donor identities. Each reaggregation group contains one existing pseudobulk; counts and contrasts remain unchanged.

## Independent reference calculation

Use the module's pinned edgePython 0.2.6 environment. Download `pbmc/combined_filtered_pca.h5ad` from the normal test-data base. The file used to establish these expectations has SHA256 `d6a8c8d837596efece4a44411d413727b8d2f568a39116ad1f260bbb2767f358`.

The reference calculation aggregates the original single-cell matrix directly with NumPy, independently of the production decoupler aggregation and the pipeline's formula selection, coefficient lookup and Parquet writer:

```python
import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd

x = ad.read_h5ad("combined_filtered_pca.h5ad")
x = x[x.obs.label == "TypeA"].copy()
samples = x.obs.groupby(["donor", "condition"], observed=True, sort=True).indices
metadata = pd.DataFrame(list(samples), columns=["donor", "condition"])
counts = np.column_stack([
    np.asarray(x.X[indices].sum(axis=0)).ravel()
    for indices in samples.values()
])
treatment = (metadata.condition == "Disease").astype(float)
paired = np.column_stack([
    np.ones(len(metadata)),
    (metadata.donor == "D02").astype(float),
    (metadata.donor == "D03").astype(float),
    treatment,
])
unpaired = np.column_stack([np.ones(len(metadata)), treatment])
for design in (paired, unpaired):
    y = ep.make_dgelist(
        counts=counts, samples=metadata,
        genes=pd.DataFrame({"genes": x.var_names}),
    )
    y = ep.calc_norm_factors(y)
    y = ep.estimate_disp(y, design=design)
    fit = ep.glm_ql_fit(y, design)
    result = ep.top_tags(
        ep.glm_ql_ftest(fit, coef=design.shape[1] - 1), n=x.n_vars,
    )["table"]
    print(result[result.genes.isin(["FOXN3", "H4C15"])])
```

| Design   | Gene  | Expected log2fc      | Expected pvalue      | Expected padj     |
| -------- | ----- | -------------------- | -------------------- | ----------------- |
| Paired   | FOXN3 | 1.2295950441713215   | 0.008705917844970962 | 1.0               |
| Unpaired | H4C15 | -0.48054392608872826 | 0.006511532640851087 | 0.928263826539495 |

With the original `estimate_disp(y)` call, the paired FOXN3 pvalue is `0.00877228752538778` and the unpaired H4C15 pvalue is `0.006495401875202621`. Both differ from the expected value by more than the `1e-7` absolute tolerance. The fold-change tolerance is `1e-6`. These checks test design consistency within the pinned library, not agreement between independently implemented statistical packages.
