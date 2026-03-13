# nf-core/scdownstream: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [2024-10-17]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `singleR` module for automated cell type annotation.
- Added `scDblFinder` module for doublet detection.
- Added optional `doublet_rate` column in input samplesheet to provide per-sample expected doublet rate for `scDblFinder`.

### `Fixed`

- Updated `scDblFinder` to use internal `dbr` estimation when `doublet_rate` is not provided, and to use provided `doublet_rate` when available.

### `Dependencies`

### `Deprecated`
