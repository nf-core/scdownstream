<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-scdownstream_logo_dark.png">
    <img alt="nf-core/scdownstream" src="docs/images/nf-core-scdownstream_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/scdownstream/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/scdownstream/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/scdownstream/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/scdownstream/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/scdownstream/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/scdownstream)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23scdownstream-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/scdownstream)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/scdownstream** is a bioinformatics pipeline that can be used to process already quantified single-cell RNA-seq data. It takes a samplesheet and h5ad-, SingleCellExperiment/Seurat- or CSV files as input and performs quality control, integration, dimensionality reduction and clustering. It produces an integrated h5ad and SingleCellExperiment file and an extensive QC report.

The pipeline is based on the learnings and implementations from the following pipelines (alphabetical):

- [panpipes](https://github.com/DendrouLab/panpipes)
- [scFlow](https://combiz.github.io/scFlow/)
- [scRAFIKI](https://github.com/Mye-InfoBank/scRAFIKI)
- [YASCP](https://github.com/wtsi-hgi/yascp)

# ![nf-core/scdownstream](docs/images/metromap.png)

Steps marked with the boat icon are not yet implemented. For the other steps, the pipeline uses the following tools:

1. Per-sample preprocessing
   1. Convert all RDS files to h5ad format
   2. Present QC for raw counts ([`MultiQC`](http://multiqc.info/))
   3. Remove ambient RNA
      - [decontX](https://bioconductor.org/packages/release/bioc/html/decontX.html)
      - [soupX](https://cran.r-project.org/web/packages/SoupX/readme/README.html)
      - [cellbender](https://cellbender.readthedocs.io/en/latest/)
      - [scAR](https://docs.scvi-tools.org/en/stable/user_guide/models/scar.html)
   4. Apply user-defined QC filters (can be defined per sample in the samplesheet)
   5. Doublet detection (Majority vote possible)
      - [SOLO](https://docs.scvi-tools.org/en/stable/user_guide/models/solo.html)
      - [scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)
      - [DoubletDetection](https://doubletdetection.readthedocs.io/en/v2.5.2/doubletdetection.doubletdetection.html)
      - [SCDS](https://bioconductor.org/packages/devel/bioc/vignettes/scds/inst/doc/scds.html)
2. Sample aggregation
   1. Merge into a single h5ad file
   2. Present QC for merged counts ([`MultiQC`](http://multiqc.info/))
   3. Integration
      - [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)
      - [scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html)
      - [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html)
      - [BBKNN](https://github.com/Teichlab/bbknn)
      - [Combat](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pp.combat.html)
      - [Seurat](https://satijalab.org/seurat/articles/integration_introduction)
3. Cell type annotation
   - [celltypist](https://www.celltypist.org/)
4. Clustering and dimensionality reduction
   1. [Leiden clustering](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)
   2. [UMAP](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html)

### Compatibility

These algorithms are sensitive at times. When sequencing depths are low, too many 0s for a particular gene or cell are likely to cause runtime errors.
That is why cells and genes can be filtered from the input that are not showing a minimal abundance. But some algorithms also explicitly demand the
unfiltered "raw" files.

The nf-core/scrnaseq pipeline offers a series of aligners, see https://nf-co.re/scrnaseq/2.7.0/parameters/ . The default --aligner "alevin" only provides raw data.
This may cause such problems with an early step in scdownstream, the removal of ambient RNA, for which the default method "decontX" expects filtered input.

The descriptions of the usage and parameters provide further guidance.
With some confidence, over time the default parameters for the nf-core scrnaseq and scdownstream will strive for maximal robustness.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

```csv title="samplesheet.csv"
sample,file
sample1,/absolute/path/to/sample1.h5ad
sample2,/absolute/path/to/sample3.h5
sample3,relative/path/to/sample2.rds
sample4,/absolute/path/to/sample3.csv
```

Each row represents a h5ad, h5, RDS or CSV file. RDS files may contain any object that can be converted to a SingleCellExperiment using the [Seurat `as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.singlecellexperiment) function.
CSV files should contain a matrix with genes as columns and cells as rows. The first column should contain cell names/barcodes.

-->

Now, you can run the pipeline using:

```bash
nextflow run nf-core/scdownstream \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/scdownstream/usage) and the [parameter documentation](https://nf-co.re/scdownstream/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/scdownstream/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/scdownstream/output).

## Credits

nf-core/scdownstream was originally written by [Nico Trummer](https://github.com/nictru).

We thank the following people for their extensive assistance in the development of this pipeline (alphabetical):

- [Fabian Rost](https://github.com/fbnrst)
- [Fabiola Curion](https://github.com/bio-la)
- [Gregor Sturm](https://github.com/grst)
- [Jonathan Talbot-Martin](https://github.com/jtalbotmartin)
- [Lukas Heumos](https://github.com/zethson)
- [Matiss Ozols](https://github.com/maxozo)
- [Nathan Skene](https://github.com/NathanSkene)
- [Nurun Fancy](https://github.com/nfancy)
- [Riley Grindle](https://github.com/Riley-Grindle)
- [Ryan Seaman](https://github.com/RPSeaman)
- [Steffen Möller](https://github.com/smoe)
- [Wojtek Sowinski](https://github.com/WojtekSowinski)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#scdownstream` channel](https://nfcore.slack.com/channels/scdownstream) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/scdownstream for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
