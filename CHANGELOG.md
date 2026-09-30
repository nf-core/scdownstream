# nf-core/scdownstream: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Changed`

- Migrate local modules to nf-core container metadata, Wave images and Conda lock files by @nictru and Codex [[#311](https://github.com/nf-core/scdownstream/pull/311)].
- Write H5AD strings in the legacy encoding pipeline-wide via `ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=0`, replacing per-module workarounds, and move `ADATA_MYGENE`, `ADATA_SETINDEX`, `ADATA_UNIFY`, `ADATA_EXTEND` and `ADATA_SPLITCOL` to anndata 0.13 [[#311](https://github.com/nf-core/scdownstream/pull/311)].
- Split nf-test CI into a module/subworkflow tier on 4-CPU runners and a pipeline tier on 16-CPU runners with one pipeline test per runner [[#317](https://github.com/nf-core/scdownstream/pull/317)].
- Split the nf-test CI module/subworkflow tier into separate `modules` and `subworkflows` categories, and run all three categories from a single `nf-test` matrix job by @nictru and Cursor [[#319](https://github.com/nf-core/scdownstream/pull/319)].
- Use `log1p` normalisation in the per-label pipeline nf-tests, which spent about 10 minutes each in `SCRAN_NORMALIZATION` on the 32k-cell extension base, by @nictru and Cursor [[#319](https://github.com/nf-core/scdownstream/pull/319)].
- Fold the `main_pipeline_build` pipeline test into the `default` test: `-profile test` now also runs Seurat integration and uses the `Adult_COVID19_PBMC` CellTypist model to match the PBMC test data [[#317](https://github.com/nf-core/scdownstream/pull/317)].
- Size CPU and memory per process for `-profile test` pipeline runs (`conf/test_resources.config`), based on measured peak usage, so that several tasks share a runner and the pipeline nf-tests run faster [[#317](https://github.com/nf-core/scdownstream/pull/317)].
- Parallelise `computeSumFactors` in `SCRAN_NORMALIZATION` over `task.cpus` and use a smaller input for its module test [[#317](https://github.com/nf-core/scdownstream/pull/317)].
- Set `KMP_AFFINITY`, `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` pipeline-wide in the `env` scope of `nextflow.config` instead of in each local Python template by @nictru and Cursor [[#318](https://github.com/nf-core/scdownstream/pull/318)].
- Move agent conventions from `AGENTS.md` to the pipeline-specific contribution guidelines in `docs/CONTRIBUTING.md` and document pipeline conventions there by @nictru and Cursor [[#318](https://github.com/nf-core/scdownstream/pull/318)].
- Replace the deprecated `QUARTONOTEBOOK` module with `QUARTO_NOTEBOOK` for rendering the QC report.
- Make final AnnData-to-RDS conversion (`ADATA_TORDS`) opt-in via `--tords` (previously always run).
- Restructure `06_per_group` to a context-first hierarchy (`{integration}/{subset}/{leiden|label}/...`) and move former `07_pseudobulk_de` outputs under `06_per_group/.../differential_expression/`.
- Prefix published result directories with numeric stage IDs (`01_load_h5ad` through `11_multiqc`) so lexical sort matches pipeline order; nest gene unification under `02_quality_control/unify`.
- Align MultiQC section order with the same pipeline sequence via `report_section_order` and `custom_content.order` in `assets/multiqc_config.yml`.
- Group raw and unified gene UpSet plots under a shared MultiQC parent section (`genes_upset`).
- Split volcano plotting out of the DE engines into a shared `custom/volcanoplot` module that reads standardised `*_results.parquet` tables [[#315](https://github.com/nf-core/scdownstream/pull/315)].
- Group differential expression volcano plots under method-specific MultiQC parents (`{integration}: {method}`) and include the DE method (including Scanpy statistical tests) in section titles [[#315](https://github.com/nf-core/scdownstream/pull/315)].
- Use contrast-first MultiQC volcano titles: Scanpy `{groupby}={group} vs rest` with optional within-filter, and PyDESeq2 / edgePython / edgepython_sc `{treatment} vs {reference} (within celltype=...)` [[#315](https://github.com/nf-core/scdownstream/pull/315)].
- Collapse Scanpy `rank_genes_groups` volcanoes into one multi-panel figure per comparison scope instead of one PNG per group [[#315](https://github.com/nf-core/scdownstream/pull/315)].
- For two-group Scanpy comparisons, show a single `{A} vs {B}` volcano instead of mirrored `{A} vs rest` and `{B} vs rest` panels [[#315](https://github.com/nf-core/scdownstream/pull/315)].
- Read `force_obs_cols` in `ADATA_MERGE` from `ext.force_obs_cols` instead of `params` by @nictru and Claude [[#XXX](https://github.com/nf-core/scdownstream/pull/XXX)].
- Pass `scvi_max_epochs` to `AMBIENT_CORRECTION` as a `take:` input instead of reading `params.scvi_max_epochs` by @nictru and Claude [[#XXX](https://github.com/nf-core/scdownstream/pull/XXX)].

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
- Deduplicate by-sample LIANA obs column selection so `context_key=condition` does not build a 2D frame and crash pandas.
- Sort doublet prediction columns before writing to `obs` so nf-test snapshots are order-stable across parallel methods.
- Reject samplesheets with duplicate `sample` names during input validation by @nictru and Claude [[#325](https://github.com/nf-core/scdownstream/pull/325)].

- Make by-sample LIANA subsampling context/donor-aware and drop contexts with too few cells per group before calling LIANA.
- Updated `scDblFinder` to use internal `dbr` estimation when `doublet_rate` is not provided, and to use provided `doublet_rate` when available.

### `Dependencies`

### `Deprecated`
