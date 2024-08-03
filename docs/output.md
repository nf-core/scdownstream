# nf-core/scdownstream: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. Per-sample preprocessing
   1. Convert all RDS files to h5ad format
   2. Create filtered matrix (if not provided)
   3. Present QC for raw counts ([`MultiQC`](http://multiqc.info/))
   4. Remove ambient RNA
      - [decontX](https://bioconductor.org/packages/release/bioc/html/decontX.html)
      - [soupX](https://cran.r-project.org/web/packages/SoupX/readme/README.html)
      - [cellbender](https://cellbender.readthedocs.io/en/latest/)
      - [scAR](https://docs.scvi-tools.org/en/stable/user_guide/models/scar.html)
   5. Apply user-defined QC filters (can be defined per sample in the samplesheet)
   6. Doublet detection (Majority vote possible)
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

### Per-sample preprocessing

<details markdown="1">
<summary>Output files</summary>

- `preprocess/${sample_id}/`
  - `converted/`: Contains the result of converting input file formats to h5ad.
  - `unified/`: Versions of the input files that have been optimized for usage in the pipeline.
  - `empty_droplet_removal/`: Results of empty droplet removal. Only if no `filtered` matrix is provided in the samplesheet.
  - `qc_raw/`: QC plots for the raw input data.
  - `ambient_rna_removal/`: Results of ambient RNA removal.
  - `custom_thresholds/`: Results of applying user-defined QC thresholds.
  - `doublet_detection/`: Directories related to doublet detection.
    - `input_rds/`: RDS version of the h5ad file that is used as input to the doublet detection tools.
    - `(doubletdetection|scds|scrublet|solo)/`: Results of doublet detection. Each directory contains a filtered `h5ad`/`rds` and a `csv`/`pkl` file with the doublet annotations.
    - `${sample_id}.h5ad`: The h5ad without doublets.
  - `qc_preprocessed/`: QC plots for the preprocessed data.

</details>

`nf-core/scdownstream` covers a range of preprocessing methods. The output of each step is stored in the `preprocess` directory. The `preprocess` directory contains a subdirectory for each sample, which contains the results of the preprocessing steps.

### Sample aggregation

<details markdown="1">
<summary>Output files</summary>

- `combine/`
  - `merge/`
    - `merged_inner.h5ad`: The merged h5ad file with only the intersection of the genes. Will be used for integration.
    - `merged_outer.h5ad`: The merged h5ad file with all genes. Will be used as base for the final h5ad file.
    - `merged_sample_genes.png`: UpSet plot showing the overlap of genes between samples.
  - `integrate/`
    - `input_hvg`
      - `*.h5ad`: The h5ad file that is used as input to the integration tools.
      - `*.rds`: RDS version of the h5ad file.
    - `${tool}`
      - `*.h5ad/*.rds`: The integrated h5ad or rds file.
      - `X_${tool}.pkl`: Low-dimensional representation of the integrated data.

</details>

The `combine` directory contains the results of the sample aggregation step. The `merge` directory contains the merged h5ad files, which are used as input to the integration tools. The `integrate` directory contains the results of the integration step. The integrated h5ad files are stored in subdirectories named after the integration tool used.

### Cell type annotation

<details markdown="1">
<summary>Output files</summary>

- `celltypes/`
  - `celltypist/`
    - `*.h5ad`: The h5ad file with cell type annotations.
    - `*.pkl`: The cell type annotations in a pickle file.

</details>

The `celltypes` directory contains the results of the cell type annotation step. So far, only `celltypist` is supported.

### Clustering and dimensionality reduction

<details markdown="1">
<summary>Output files</summary>

- `cluster_dimred/`
  - `${integration}/`
    - `neighbors/`
      - `*.h5ad`: The h5ad file with the neighborhood graph.
    - `leiden/`
      - ${resolution}/`
        - `*.h5ad`: The h5ad file with the leiden clustering.
        - `*.pkl`: The leiden clustering in a pickle file.
    - `umap/`
      - `*.h5ad`: The h5ad file with the UMAP coordinates.
      - `*.pkl`: The UMAP coordinates in a pickle file.

</details>

The `cluster_dimred` directory contains the results of the clustering and dimensionality reduction step. The results are stored in subdirectories named after the integration tool used.

### Finalize

<details markdown="1">
<summary>Output files</summary>

- `finalize/`
  - `merged.h5ad`: The final h5ad file with all results.
  - `merged.rds`: RDS version of the final h5ad file.
  - `merged_metadata.csv`: Metadata of the final h5ad file.

</details>

The `finalize` directory contains the final results of the pipeline. The final h5ad file contains all results from the pipeline and is stored in the `merged.h5ad` file. The metadata of the final h5ad file is stored in the `merged_metadata.csv` file.

### MultiQC

:::warning
The MultiQC report is not yet fully implemented and will be improved in future releases.
Until then, feel free to utilize [CELLxGENE](https://cellxgene.cziscience.com/docs/05__Annotate%20and%20Analyze%20Your%20Data/5_1__Getting%20Started:%20Install,%20Launch,%20Quick%20Start) for interactive exploration of the results.
:::

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
