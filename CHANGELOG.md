# nf-core/scdownstream: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [2024-10-17]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Changed`

- Make final AnnData-to-RDS conversion (`ADATA_TORDS`) opt-in via `--tords` (previously always run).
- Restructure `06_per_group` to a context-first hierarchy (`{integration}/{subset}/{leiden|label}/...`) and move former `07_pseudobulk_de` outputs under `06_per_group/.../differential_expression/`.
- Prefix published result directories with numeric stage IDs (`01_load_h5ad` through `11_multiqc`) so lexical sort matches pipeline order; nest gene unification under `02_quality_control/unify`.
- Align MultiQC section order with the same pipeline sequence via `report_section_order` and `custom_content.order` in `assets/multiqc_config.yml`.
- Group raw and unified gene UpSet plots under a shared MultiQC parent section (`genes_upset`).
- Group differential expression volcano plots under method-specific MultiQC parents (`{integration}: {method}`) and include the DE method (including Scanpy statistical tests) in section titles.
- Use contrast-first MultiQC volcano titles: Scanpy `{groupby}={group} vs rest` with optional within-filter, and PyDESeq2 / edgePython / edgepython_sc `{treatment} vs {reference} (within celltype=...)`.
- Collapse Scanpy `rank_genes_groups` volcanoes into one multi-panel figure per comparison scope instead of one PNG per group.
- For two-group Scanpy comparisons, show a single `{A} vs {B}` volcano instead of mirrored `{A} vs rest` and `{B} vs rest` panels.

### `Added`

- Add opt-in `--tords` to convert the final AnnData object to RDS via `ADATA_TORDS` (off by default).
- Add LIANA rank-aggregate dotplot, circle, and tileplot PNGs with MultiQC embedding.
- Add opt-in Tensor-cell2cell analysis: by-sample LIANA (`liana/bysample`) followed by tensor factorisation and sender-receiver loadings-product heatmaps (`cell2cell/tensor`).
- Add volcano plots and optional `--interesting_genes` highlighting for Scanpy, PyDESeq2, edgePython, and edgepython_sc differential expression.
- Add CyteType module for automated cell type annotation.
- Add ribosomal/haemoglobin QC metrics [[#277]https://github.com/nf-core/scdownstream/pull/277]
- Add reporting using Quarto [[#258](https://github.com/nf-core/scdownstream/pull/258)]
- Convert to Nextflow strict mode [[#244](https://github.com/nf-core/scdownstream/pull/244)]
- Add `singleR` module for automated cell type annotation [[#200](https://github.com/nf-core/scdownstream/pull/200)]
- Use topics for software versioning [[#252](https://github.com/nf-core/scdownstream/pull/252)]
- Added `singleR` module for automated cell type annotation.
- Added `scDblFinder` module for doublet detection.
- Added optional `doublet_rate` column in input samplesheet to provide per-sample expected doublet rate for `scDblFinder`.

### `Fixed`

- Pin `pandas` and `anndata` in `SCANPY_RANKGENESGROUPS`, `CYTETYPE`, and `ADATA_EXTEND` to avoid nullable-string H5AD writes; remove `allow_write_nullable_strings` opt-ins and cast remaining `StringDtype` columns to plain `object` at write boundaries.
- Pass gene symbols into edgePython DE so volcano labels use gene names instead of row indices.

- Make by-sample LIANA subsampling context/donor-aware and drop contexts with too few cells per group before calling LIANA.
- Updated `scDblFinder` to use internal `dbr` estimation when `doublet_rate` is not provided, and to use provided `doublet_rate` when available.

### `Dependencies`

### `Deprecated`
