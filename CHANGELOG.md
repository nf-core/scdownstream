# nf-core/scdownstream: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [2024-10-17]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Changed`

- Prefix published result directories with numeric stage IDs (`01_load_h5ad` through `11_multiqc`) so lexical sort matches pipeline order; nest gene unification under `02_quality_control/unify`.
- Align MultiQC section order with the same pipeline sequence via `report_section_order` and `custom_content.order` in `assets/multiqc_config.yml`.
- Group raw and unified gene UpSet plots under a shared MultiQC parent section (`genes_upset`).

### `Added`

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

- Updated `scDblFinder` to use internal `dbr` estimation when `doublet_rate` is not provided, and to use provided `doublet_rate` when available.

### `Dependencies`

### `Deprecated`
