# Single-cell fold-change regression

`units.nf.test` executes the real module on the existing combined fixture. `TypeA` contains 384 cells, three donors and two conditions. The reference condition is `Healthy` and the tested coefficient is `Disease`.

## Independent reference calculation

Use the module's pinned edgePython 0.2.6 environment. The SHA256 of `pbmc/combined_filtered_pca.h5ad` used for this calculation is `d6a8c8d837596efece4a44411d413727b8d2f568a39116ad1f260bbb2767f358`.

Run the library directly, without the process template, gene-label insertion or standard-format writer:

```python
import anndata as ad
import edgepython as ep
import numpy as np

x = ad.read_h5ad("combined_filtered_pca.h5ad")
x = x[x.obs.label == "TypeA"].copy()
design = np.column_stack([
    np.ones(x.n_obs), (x.obs.condition == "Disease").astype(float),
])
fit = ep.glm_sc_fit(
    x.X.toarray().T, design=design,
    sample=x.obs.donor.astype(str).to_numpy(), norm_method="TMM",
)
fit = ep.shrink_sc_disp(fit, robust=True)
genes = np.asarray(x.var_names)[fit["gene_mask"]]
result = ep.glm_sc_test(fit, coef=1)["table"]
for gene in ("FOXN3", "CRTAP"):
    i = np.flatnonzero(genes == gene)[0]
    print(gene, fit["coefficients"][i, 1] / np.log(2), result.iloc[i])
```

The model uses a natural-log link. A treatment coefficient `beta` represents an expression ratio `exp(beta)`; its base-2 logarithm is `beta / log(2)`. In pinned edgePython, `glm_sc_test` places `fit['coefficients'][:, coef]` directly into `logFC`, without conversion. The expected base-2 coefficients are `1.2370226089406087` for FOXN3 and `-2.062184319422321` for CRTAP. The original process writes `0.8574387336760907` and `-1.4293972468025113`, respectively. The `1e-4` absolute tolerance distinguishes both the original units and double conversion.

The regression also checks the unchanged FOXN3 pvalue (`0.004165649756740918`, tolerance `1e-6`) and adjusted pvalue (`0.9614059396460176`, tolerance `1e-5`), gene count, schema and contrast metadata. It validates the writer against the pinned model coefficients; it does not independently validate the model's statistical estimator.
