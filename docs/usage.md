# nf-core/scdownstream: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/scdownstream/usage](https://nf-co.re/scdownstream/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Filtered and unfiltered matrices

Throughout this documentation, you will find references to `filtered` and `unfiltered` matrices.
The `unfiltered` matrices are matrices which still contain empty droplets, whereas the `filtered` matrices have been filtered for empty droplets.
A more technical definition can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices).
`CellRanger` provides you with both matrices, whereas other quantification tools only provide you with the `unfiltered` matrix.
The pipeline can handle the following cases:

1. You have both `filtered` and `unfiltered` matrices: Provide both matrices in the samplesheet and the pipeline will use the `unfiltered` matrix for ambient RNA removal and the `filtered` matrix for all other steps.
2. You only have the `filtered` matrix: Provide the `filtered` matrix in the samplesheet and the pipeline will use it for all steps.
   SoupX is the default ambient correction method and requires an unfiltered matrix. For filtered-only input, disable ambient correction per sample (`ambient_correction=false`) or set `--ambient_correction decontx`.
3. You only have the `unfiltered` matrix: Provide the `unfiltered` matrix in the samplesheet and the pipeline will automatically create a `filtered` matrix by identifying empty droplets using `CellBender`.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
Use this parameter to specify its location.
It has to be a comma-separated file with at least 2 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Minimal samplesheet

The samplesheet needs to contain at least two columns: `sample` and at least one out of `filtered` and `unfiltered`:

```csv title="samplesheet.csv"
sample,unfiltered
sample1,/absolute/path/to/sample1.h5ad
sample2,relative/path/to/sample2.rds
sample3,/absolute/path/to/sample3.csv
```

### Full samplesheet

There are a couple of optional columns that can be used for more advanced features:

```csv title="samplesheet.csv"
sample,filtered,unfiltered,batch_col,label_col,condition_col,unknown_label,min_genes,min_cells,min_counts_cell,min_counts_gene,expected_cells,doublet_rate,ambient_correction,ambient_corrected_integration
sample1,/absolute/path/to/sample1_filtered.h5ad,/absolute/path/to/sample1.h5ad,batch,cell_type,condition,unknown,1,2,3,4,5000,0.08,true,false
sample2,relative/path/to/sample2_filtered.rds,relative/path/to/sample2.rds,batch_id,annotation,condition,unannotated,5,6,7,8,3000,,false,
sample3,/absolute/path/to/sample3_filtered.csv,/absolute/path/to/sample3.csv,,,,,9,10,11,12,,,true,true
```

For CSV input files, specifying the `batch_col`, `label_col`, `condition_col`, and `unknown_label` columns will not have any effect, as no additional metadata is available in the CSV file.

| Column                             | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                           | Unique sample identifier. Will be added to the pipeline output objects as `sample` column.                                                                                                                                                                                                                                                                                                                                                                 |
| `filtered`                         | May contain paths to `h5ad`, `h5`, `rds`, or `csv` files. `rds` files may contain any object that can be converted to a `SingleCellExperiment` using the [Seurat `as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.singlecellexperiment) function. `csv` files should contain a matrix with genes as columns and cells as rows.                                                                                                         |
| `unfiltered`                       | Same as `filtered`, but for the unfiltered cellranger or nf-core/scrnaseq output. If not provided, only `decontX` can be used for ambient RNA removal.                                                                                                                                                                                                                                                                                                     |
| `batch_col`                        | Column in the input file containing batch information. If not provided, the entire input object will be considered as one batch. If the `batch_col` is something else than `batch`, it will be renamed to `batch` during pipeline execution.                                                                                                                                                                                                               |
| `symbol_col`                       | Column in the input file containing gene symbol information. Defaults to `index`. There are two special values that can be used: `index` and `none`. `index` will use the row names of the matrix as gene symbols. `none` will trigger the pipeline to perform gene symbol conversion using MyGene.info based on the `geneid_col` and the pipeline `--species` parameter. The values from `symbol_col` will be set as the index during pipeline execution. |
| `geneid_col`                       | Column in the input file containing gene identifier information. Defaults to `index`. Only used if `symbol_col` is set to `none`.                                                                                                                                                                                                                                                                                                                          |
| `label_col`                        | Column in the input file containing cell type information. Defaults to `label`. If the column does not exist in the input object, the pipeline will create a new column and put `unknown` in it. If the `label_col` is something else than `label`, it will be renamed to `label` during pipeline execution.                                                                                                                                               |
| `condition_col`                    | Column in the input file containing condition information (e.g. disease state, treatment). If the column does not exist in the input object, the pipeline will create a new column and put `unknown` in it. If the `condition_col` is something else than `condition`, it will be renamed to `condition` during pipeline execution.                                                                                                                        |
| `donor_col`                        | Column in the input file containing biological replicate / donor identifiers (e.g. patient, mouse). Required in the samplesheet when pseudobulking is enabled. If the column is something else than `donor`, it will be renamed to `donor` during unification.                                                                                                                         |
| `unknown_label`                    | Value in the `label_col` column that should be considered as unknown. Defaults to `unknown`. If the `unknown_label` is something else than `unknown`, it will be renamed to `unknown` during pipeline execution. If trying to perform integration with scANVI, more than one unique label other than `unknown` must exist in the input data.                                                                                                               |
| `counts_layer`                     | Layer in the input file containing the raw counts matrix. Defaults to `X`.                                                                                                                                                                                                                                                                                                                                                                                 |
| `min_genes`                        | Minimum number of genes required for a cell to be considered. Defaults to `0`.                                                                                                                                                                                                                                                                                                                                                                             |
| `min_cells`                        | Minimum number of cells required for a gene to be considered. Defaults to `20`.                                                                                                                                                                                                                                                                                                                                                                             |
| `min_counts_cell`                  | Minimum number of counts required for a cell to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                                                            |
| `min_counts_gene`                  | Minimum number of counts required for a gene to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                                                            |
| `expected_cells`                   | Number of expected cells, used as input to CellBender for empty droplet detection.                                                                                                                                                                                                                                                                                                                                                                         |
| `doublet_rate`                     | Optional expected doublet rate (0-1) for `scDblFinder`. If not provided, `scDblFinder` estimates it internally.                                                                                                                                                                                                                                                                                                                                            |
| `max_mito_percentage`              | Maximum percentage of mitochondrial reads for a cell to be considered. Defaults to `8`.                                                                                                                                                                                                                                                                                                                                                                  |
| `min_ribo_percentage`              | Minimum percentage of ribosomal reads for a cell to be considered. Defaults to `0`.                                                                                                                                                                                                                                                                                                                                                                        |
| `max_hb_percentage`                | Maximum percentage of haemoglobin reads for a cell to be considered. Defaults to `100`.                                                                                                                                                                                                                                                                                                                                                                    |
| `log1p_total_counts_nmads`         | MAD cutoff for `log1p_total_counts`. Cells outside `median ± n × MAD` are removed. Defaults to `5`. Set to `0` to disable.                                                                                                                                                                                                                                                                                                                                             |
| `log1p_n_genes_by_counts_nmads`    | MAD cutoff for `log1p_n_genes_by_counts`. Defaults to `5`. Set to `0` to disable.                                                                                                                                                                                                                                                                                                                                                                                      |
| `pct_counts_in_top_20_genes_nmads` | MAD cutoff for `pct_counts_in_top_20_genes`. Defaults to `5`. Set to `0` to disable.                                                                                                                                                                                                                                                                                                                                                                                   |
| `pct_counts_mt_nmads`              | MAD cutoff for `pct_counts_mt`. Defaults to `3`. Set to `0` to disable.                                                                                                                                                                                                                                                                                                                                                                                                |
| `ambient_correction`               | Whether to perform ambient RNA correction for this sample. Set to `true` to use the globally configured method, `false` to skip ambient correction for this sample. Defaults to `true`.                                                                                                                                                                                                                                                                    |
| `ambient_corrected_integration`    | Whether to use ambient-corrected counts for integration for this sample. Set to `true` to use corrected counts in downstream integration, `false` to store them only as additional layers. Can override the global `--ambient_corrected_integration` parameter. Defaults to global setting.                                                                                                                                                                |

MAD-based filtering follows the [sc-best-practices recipe](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-cells). QC thresholds are defined in the [input schema](../assets/schema_input.json) and applied per sample when the samplesheet is read; omitted columns receive those schema defaults. Each MAD column is independent; set a value to `0` to disable that metric for a sample.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### sc-best-practices defaults

Pipeline defaults match the [sc-best-practices QC chapter](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html) out of the box. A minimal samplesheet with only `sample` and matrix paths is sufficient:

```csv title="samplesheet.csv"
sample,filtered
sample1,/path/to/sample1_filtered.h5ad
```

With no QC columns, each sample receives:

- MAD filtering: `log1p_total_counts_nmads=5`, `log1p_n_genes_by_counts_nmads=5`, `pct_counts_in_top_20_genes_nmads=5`, `pct_counts_mt_nmads=3`
- `max_mito_percentage=8`
- `min_cells=20` (gene filter after ambient correction)

Doublet handling also follows the book: `--doublet_detection` defaults to `scdblfinder`, and `--doublet_removal` defaults to `false`, so doublet scores are written to `adata.obs` without removing cells until you opt in. Ambient RNA correction defaults to SoupX and requires an unfiltered matrix for each corrected sample.

Before integration, merged raw counts are subset to informative genes. By default, [`feature_selection`](https://nf-co.re/scdownstream/parameters#feature_selection) is `deviance` (binomial deviance via scry, ~4,000 genes when [`integration_n_features`](https://nf-co.re/scdownstream/parameters#integration_n_features) is `0`). Use `--feature_selection hvgs` for scanpy highly variable genes, or `--feature_selection none` to skip gene filtering. The `python_only` profile sets `feature_selection` to `hvgs` automatically.

Optional normalisation layers can be computed before feature selection with [`normalization_methods`](https://nf-co.re/scdownstream/parameters#normalization_methods) (`log1p`, `scran`, `pearson_residuals`). Raw counts remain in `X`; shifted-log output is stored in `log1p_norm`, scran output in `scran`, and analytic Pearson residuals in `pearson_residuals`. Use [`transformed_layer`](https://nf-co.re/scdownstream/parameters#transformed_layer) to point PCA integration (and HVG selection when applicable) at one of these layers. Count-model integrations (`scvi`, `scanvi`) and all differential expression engines continue to use raw counts.

To disable individual filters for a sample, set the corresponding samplesheet column to `0` (for MAD metrics or `min_cells`) or a permissive value (e.g. `max_mito_percentage=100`):

```csv title="samplesheet.csv"
sample,filtered,log1p_total_counts_nmads,pct_counts_mt_nmads,max_mito_percentage,min_cells
sample1,/path/to/sample1.h5ad,0,0,100,0
```

To remove detected doublets after inspection, pass `--doublet_removal true` (optionally with `--doublet_detection_threshold` to require agreement across multiple tools).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scdownstream --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile.
See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors.
> Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scdownstream -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Cell type annotation

The pipeline supports automated cell type annotation with [Celltypist](https://github.com/Teichlab/celltypist), [singleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html), and [CyteType](https://github.com/NygenAnalytics/cytetype). Each method is optional and controlled by its own parameters below.

#### Celltypist

[Celltypist](https://github.com/Teichlab/celltypist) annotates cells using pretrained logistic regression models. Specify the models with [`celltypist_model`](https://nf-co.re/scdownstream/dev/parameters/#celltypist_model); use a comma-separated list for multiple models. Available models are listed on the [Celltypist models page](https://www.celltypist.org/models). When this parameter is empty (the default), Celltypist is skipped.

#### singleR

For `singleR`, you can provide a CSV file with information about the celldex references to use for the singleR cell type annotation with the [`celldex_reference` parameter](https://nf-co.re/scdownstream/dev/parameters/#celldex_reference).
The exising references are described in the [celldex package description](https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf).
You can also provide paths to tar archives of pre-downloaded references (useful if your runtime environment does not have access to the internet).

A CSV file that refers to the celldex references via name can look like this:

```csv title="celldex_references.csv"
id,label,reference,version
hpca,label.main,hpca,2024-02-26
monaco_immune,label.fine,monaco_immune,2024-02-26
```

A CSV file that refers to the celldex references via path can look like this:

```csv title="celldex_references.csv"
hpca,label.main,/path/to/hpca.tar
monaco_immune,label.fine,/path/to/monaco_immune.tar
```

Example tar archives can be found [here](https://github.com/nf-core/test-datasets/tree/scdownstream/singleR).

#### CyteType

[CyteType](https://github.com/NygenAnalytics/cytetype) is a multi-agent LLM-driven annotator that takes per-cluster marker genes and a free-text study description and returns predicted cell type labels. The pipeline runs CyteType on merged data after integration, clustering, and global differential expression — once per grouping (each Leiden resolution and label column). Cluster labels and marker genes are taken automatically from each grouping's obs column and `uns['rank_genes_groups']`.

To enable CyteType, set [`cytetype_study_context`](https://nf-co.re/scdownstream/dev/parameters/#cytetype_study_context) to a short free-text description of your study (the more specific, the better). When this parameter is empty (the default), CyteType is skipped. In the analysis plan, CyteType is controlled by the `cytetype` token. CyteType always reads Wilcoxon `rank_genes_groups` results; when `cytetype` is active, `wilcoxon` is added to the resolved `de_methods` for each eligible clustering if not already present.

```bash
nextflow run nf-core/scdownstream \
    --input samplesheet.csv \
    --outdir results \
    --cytetype_study_context "Human PBMC from healthy donor, 10X Genomics 3' scRNA-seq"
```

> [!IMPORTANT]
> CyteType calls the remote `https://cytetype.nygen.io` API and therefore **requires internet access** from the compute node running the `CYTETYPE` task.

If your CyteType deployment requires authentication, set the Nextflow secret `CYTETYPE_API_KEY` before the run (for example `nextflow secrets set CYTETYPE_API_KEY '<token>'`).

### Cell cycle scoring

Cell cycle scoring assigns each cell an S-phase score, G2M-phase score, and a predicted cell cycle phase (`S`, `G2M`, or `G1`) based on the expression of curated marker genes (Tirosh et al. 2015, same gene sets as Seurat).
The scores are stored in `adata.obs` as `S_score`, `G2M_score`, and `phase`, and are available as covariates in downstream integration steps.

Cell cycle scoring is enabled by default.
To skip it:

```bash
nextflow run nf-core/scdownstream --input samplesheet.csv --outdir results --cell_cycle_scoring false
```

#### Species

Bundled gene lists are provided for human and mouse.
`--species` also selects the MyGene.info taxonomy used when samples have `symbol_col: none` and gene identifiers are converted via MyGene.info.
Select the appropriate species with `--species`:

```bash
# mouse
nextflow run nf-core/scdownstream --input samplesheet.csv --outdir results --species mouse
```

#### Custom gene lists

For other organisms (e.g. rat, zebrafish), you can provide your own gene lists — one gene symbol per line — via `--s_genes` and `--g2m_genes`:

```bash
nextflow run nf-core/scdownstream --input samplesheet.csv --outdir results \
    --s_genes /path/to/my_s_genes.txt \
    --g2m_genes /path/to/my_g2m_genes.txt
```

The bundled gene lists can be found in [`assets/cell_cycle_genes/`](../assets/cell_cycle_genes/) and serve as templates for custom lists.

#### Using scores in downstream analysis

The `S_score` and `G2M_score` columns can be passed to integration tools as continuous covariates to regress out cell cycle effects:

```bash
nextflow run nf-core/scdownstream --input samplesheet.csv --outdir results \
    --scvi_continuous_covariates S_score,G2M_score
```

### Reference mapping and extension

**Reference mapping** means **mapping new cells into a latent space using a pre-trained model** instead of training that integration step only on the query data.
In this pipeline this can be done using **scVI**, **scANVI**, **scimilarity**, and **Symphony**.
To enable it, add the corresponding method to [`integration_methods`](https://nf-co.re/scdownstream/parameters#integration_methods) (`scvi`, `scanvi`, `scimilarity`, and/or `symphony`) and set the matching model parameters for each method you use: [`scvi_model`](https://nf-co.re/scdownstream/parameters#scvi_model), [`scanvi_model`](https://nf-co.re/scdownstream/parameters#scanvi_model), [`scimilarity_model`](https://nf-co.re/scdownstream/parameters#scimilarity_model), and [`symphony_reference`](https://nf-co.re/scdownstream/parameters#symphony_reference) (see the [parameter reference](https://nf-co.re/scdownstream/parameters) for file types, defaults, and help text).

For Symphony reference mapping, provide the compact Symphony reference AnnData from a prior de novo run (`{outdir}/combine/integrate/symphony/symphony_reference.h5ad`). It contains the gene statistics, PCA loadings, Harmony centroids, and normalization metadata required for query mapping.

**Extension** is for users that have outputs of a previous run of `nf-core/scdownstream` and want to extend it with new data, without re-running the integration from scratch.
It only works if `scvi`, `scanvi`, `scimilarity`, and/or `symphony` have been enabled in `integration_methods` in the original pipeline run.
Other integration methods than the four mentioned before are not supported for this.
In simple terms, in this setup the workflow is: (1) project new data into the latent space learned from the data in the original run, and then (2) combine the datasets.
For (1), provide the same checkpoints as for reference mapping ([`scvi_model`](https://nf-co.re/scdownstream/parameters#scvi_model), [`scanvi_model`](https://nf-co.re/scdownstream/parameters#scanvi_model), [`scimilarity_model`](https://nf-co.re/scdownstream/parameters#scimilarity_model), [`symphony_reference`](https://nf-co.re/scdownstream/parameters#symphony_reference)).
For (2), pass the integrated `.h5ad` from the original run as [`base_adata`](https://nf-co.re/scdownstream/parameters#base_adata).

Pre-trained scVI models are also shared on [scvi-hub](https://huggingface.co/scvi-tools).

### Integration benchmarking (scib-metrics)

You can run [scib-metrics](https://scib-metrics.readthedocs.io/) on each integration output by setting [`scib`](https://nf-co.re/scdownstream/parameters#scib) to `true` (default is `false`).
The step is **not** run when only [`base_adata`](https://nf-co.re/scdownstream/parameters#base_adata) and [`base_embeddings`](https://nf-co.re/scdownstream/parameters#base_embeddings) are provided without `--input`.

Metrics tables are published under `combine/integrate/scib_metrics/<method>/`, and a summary table is included in the MultiQC report.
Values are not numerically comparable to the original scIB reference implementation (see the scib-metrics documentation).
Rare batches or uninformative labels can make scores such as kBET unstable.

## Clustering and beyond

After integration, the pipeline builds a neighbour graph and UMAP for every integration output, then performs Leiden clustering and a suite of downstream analyses on each clustering result.

### Clustering

For each integration method, the pipeline:

1. Computes a **KNN neighbour graph** (using the reduced embedding, e.g. PCA or scVI latent space).
2. Generates a **UMAP** layout.
3. Runs **Leiden clustering** at every resolution listed in [`clustering_resolutions`](https://nf-co.re/scdownstream/parameters#clustering_resolutions) (default `0.25,0.5,1.0`).

Steps 1 and 2 always run for every integration and every subset (global and per-label). Step 3 is controlled by the analysis plan (see below).

**Per-label sub-clustering** — when [`cluster_per_label`](https://nf-co.re/scdownstream/parameters#cluster_per_label) is `true`, the pipeline splits the integrated object by the label column and builds a separate neighbour graph, UMAP, and Leiden clustering for each label value, in addition to the global clustering.

### Downstream analyses

For each Leiden clustering result the pipeline runs a configurable set of downstream analyses:

| Analysis          | What it does                                                                 | Skip parameter                                                                                       |
| ----------------- | ---------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| **PAGA**          | Trajectory / connectivity graph between clusters                             | —                                                                                                    |
| **LIANA**         | Ligand–receptor interaction analysis                                         | [`skip_liana`](https://nf-co.re/scdownstream/parameters#skip_liana)                                  |
| **DE**                       | Cell-level and sample-level differential expression via `de_methods`         | omit `de` from the analysis plan and/or set [`de_methods`](https://nf-co.re/scdownstream/parameters#de_methods) to an empty string |
| **Aggregate per-cell annotation** | Majority vote of per-cell SingleR/CellTypist labels per cluster (columns derived from annotator manifests) | omit `aggregate_per_cell_annotation` from the analysis plan and/or do not run per-cell annotators |
| **CyteType**                      | LLM-based cluster cell type annotation                                       | omit `cytetype` from the analysis plan and/or leave [`cytetype_study_context`](https://nf-co.re/scdownstream/parameters#cytetype_study_context) empty |

By default (no `--analysis_plan`), `paga`, `liana`, `de`, `aggregate_per_cell_annotation`, and `cytetype` run for every clustering result, subject to `skip_liana`, `de_methods`, and `cytetype_study_context` above. Sample-level pseudobulk DE runs automatically when `pydeseq2` or `edgepython` are included in the resolved `de_methods` for a clustering.

Per-cell SingleR and CellTypist annotation runs earlier in the pipeline (before sample merge) when [`celldex_reference`](https://nf-co.re/scdownstream/parameters#celldex_reference) and/or [`celltypist_model`](https://nf-co.re/scdownstream/parameters#celltypist_model) are set. `aggregate_per_cell_annotation` summarises those per-cell predictions per cluster using columns declared in each annotator's manifest. `cytetype` is an independent cluster-level annotator.

### Differential expression methods

The [`de_methods`](https://nf-co.re/scdownstream/parameters#de_methods) parameter selects which engines run when the corresponding analysis token is active:

| Method | Purpose |
| ------ | ------- |
| `wilcoxon` | Scanpy `rank_genes_groups` with the Wilcoxon test (default) |
| `t-test` | Scanpy `rank_genes_groups` with Student's t-test |
| `t-test_overestim_var` | Scanpy `rank_genes_groups` with t-test (overestimated variance) |
| `logreg` | Scanpy `rank_genes_groups` with logistic regression |
| `edgepython_sc` | Donor-aware single-cell mixed model ([edgePython](https://github.com/pachterlab/edgePython)) |
| `pydeseq2` | Sample-level PyDESeq2 on pseudobulk counts (triggers pseudobulk aggregation) |
| `edgepython` | Sample-level edgeR-style QL F-test on pseudobulk counts (triggers pseudobulk aggregation) |

Use multiple comma-separated values to compare methods in one run, for example `--de_methods wilcoxon,pydeseq2,edgepython`. To disable all DE engines globally, pass an empty value: `--de_methods ''`. DE still runs for a clustering when a matching `--analysis_plan` row supplies `de_methods`.

**Rank-genes-groups comparisons:** each selected Scanpy method (`wilcoxon`, `t-test`, `t-test_overestim_var`, `logreg`) runs `rank_genes_groups` for global cluster markers, per-condition cluster markers, and per-cluster condition contrasts when those columns are present. Multiple Scanpy methods produce parallel outputs (separate `uns` keys and plots per method). For donor-aware modelling across biological replicates, add `pydeseq2` and/or `edgepython` to `de_methods`, or use `edgepython_sc`.

**Pseudobulk settings:** aggregation follows the [sc-best-practices DGE tutorial](https://www.sc-best-practices.org/conditions/differential_gene_expression.html) using [decoupler](https://decoupler.readthedocs.io/) (`pp.pseudobulk` + `filter_samples`). Biological replicate identifiers are unified to a `donor` column during QC unification. Set the source column per sample with `donor_col` in the samplesheet (required when pseudobulking is enabled). When extending a previous run, `base_adata` must already contain a `donor` column. Minimum cells per pseudobulk sample ([`pseudobulk_min_num_cells`](https://nf-co.re/scdownstream/parameters#pseudobulk_min_num_cells), default `10`) and minimum total counts ([`pseudobulk_min_total_counts`](https://nf-co.re/scdownstream/parameters#pseudobulk_min_total_counts), default `1000`) filter low-coverage pseudobulk profiles. Sample-level DE uses a fixed `~ donor + condition` design per cell-type stratum. Set [`reference_condition`](https://nf-co.re/scdownstream/parameters#reference_condition) to choose the baseline condition for pseudobulk and `edgepython_sc` contrasts (defaults to the first level alphabetically). Use [`pseudobulk`](https://nf-co.re/scdownstream/parameters#pseudobulk) to export pseudobulk count matrices without running pseudobulk DE methods.

### Analysis plan

With many integration methods and resolutions the full downstream suite can generate a large number of tasks. The optional [`analysis_plan`](https://nf-co.re/scdownstream/parameters#analysis_plan) parameter accepts a CSV that controls exactly which Leiden resolutions are computed and which analyses run for each clustering result.

Each row in the CSV selects a subset of clusterings. **All columns are optional** — an empty cell acts as a wildcard that matches everything:

| Column        | Empty means                                                                              |
| ------------- | ---------------------------------------------------------------------------------------- |
| `integration` | match all integration methods                                                            |
| `subset`      | match all subsets (`global` and per-label)                                               |
| `resolution`  | match all resolutions (still bounded by `--clustering_resolutions`)                    |
| `analyses`    | run `paga`, `liana`, `de`, `aggregate_per_cell_annotation`, and `cytetype`         |
| `de_methods`  | use the global [`de_methods`](https://nf-co.re/scdownstream/parameters#de_methods) default |

When multiple rows match a clustering result, their `analyses` lists are **combined** (duplicates removed). If any matching row leaves `analyses` empty, all analyses run for that clustering. Clusterings that match **no** row are excluded from Leiden and all downstream analyses — but their UMAP and neighbour graph are still computed.

Example plan: full analysis on Symphony at resolution 0.5, DE-only at resolution 1.0 for every integration, and DE-only for scVI at any resolution:

```csv title="analysis_plan.csv"
integration,subset,resolution,analyses,de_methods
scvi,global,0.5,"de","wilcoxon,pydeseq2,edgepython"
,,,de,wilcoxon
```

```bash
nextflow run nf-core/scdownstream \
    --input samplesheet.csv \
    --outdir results \
    --analysis_plan analysis_plan.csv
```

:::note
Label-column analyses (PAGA / LIANA / DE run on the merged `label` column rather than on Leiden clusters) are not controlled by the analysis plan; they always run subject to the global skip flags.
:::

### Skipping integration

:::tip
This can be useful if you have assigned cell type annotations to the integrated object and want to perform further analysis based on these annotations.
:::

If you want to run tasks after the integration step without performing integration, you can provide a previous result of the pipeline as [`base_adata`](https://nf-co.re/scdownstream/parameters#base_adata).
You do not need to provide a samplesheet via the [`input`](https://nf-co.re/scdownstream/parameters#input) parameter in this case.
You also need either [`base_embeddings`](https://nf-co.re/scdownstream/parameters#base_embeddings) to reuse existing embeddings, or [`integrate_per_label`](https://nf-co.re/scdownstream/parameters#integrate_per_label) to compute new integrations independently for each group in [`base_label_col`](https://nf-co.re/scdownstream/parameters#base_label_col).
Set [`base_condition_col`](https://nf-co.re/scdownstream/parameters#base_condition_col) if your condition column is not named `condition`.

The pipeline will then re-execute the tasks after the integration step without performing integration again.
Most interestingly, the pipeline will generate cell type specific UMAPs, clusterings, and PAGA graphs, if [`clustering_per_label`](https://nf-co.re/scdownstream/parameters#clustering_per_label) is set to `true`.

If [`integrate_per_label`](https://nf-co.re/scdownstream/parameters#integrate_per_label) is enabled, [`base_label_col`](https://nf-co.re/scdownstream/parameters#base_label_col) is the split/grouping column, not necessarily the supervised cell-type label used by integration methods.
Use [`integrate_per_label_whitelist`](https://nf-co.re/scdownstream/parameters#integrate_per_label_whitelist) to restrict per-label integration to a subset of groups (comma-separated values from `base_label_col`). When omitted, integration runs for every group. Whitelist values must use the same filesystem-safe subset names as in `analysis_plan` (spaces replaced with underscores).
Because [`base_adata`](https://nf-co.re/scdownstream/parameters#base_adata) is expected to be a previous pipeline output, batch-aware methods use the standard `batch` column, and scANVI uses the standard `label` column with `Unknown` as the unlabeled category.
Subset names in `analysis_plan` match the filesystem-safe keys produced by splitting the AnnData object; spaces in group values are replaced with underscores.
Per-label integrations are treated as already split for clustering, so the pipeline creates subset-specific embedding keys such as `X_pca-SRR28679756_pca` and UMAP keys such as `X_pca-SRR28679756_umap` in the finalized base AnnData.

### GPU acceleration

:::warning{title="Experimental feature"}
This is an experimental feature and may produce errors.
If you encounter any issues, please report them on the [nf-core/scdownstream GitHub repository](https://github.com/nf-core/scdownstream/issues/new?assignees=&labels=bug&projects=&template=bug_report.yml).
:::

:::info{title="Prerequisites"}

- GPU acceleration has only been tested with Docker, Singularity and Apptainer.
  - Other container technologies might work, but have not been tested.
  - Conda is not supported.
- CUDA 12.0 or later is required.
- The NVIDIA GPUs must have a [Compute Capability](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities) of 7.0 or higher.
- ROCM is currently unsupported

:::

Tools with implemented support for GPU acceleration are:

- cellbender
- scvi-tools
  - scVI/scANVI
  - scAR
  - solo

To utilize GPU acceleration, you need to specify the `gpu` profile.
This will make the tool steps use cuda-enabled environments and it will tell the tools to use the GPU.
All processes which support GPU acceleration are marked with the `process_gpu` label.

You also need to make sure that the tasks are run on a machine with a GPU.
If all tasks are run on a machine with a GPU, no further action is needed.
If you are running the pipeline on a slurm cluster, where there is dedicated queue for GPU jobs, you need additional configuration that might look like this:

```bash
process {
  withLabel:process_gpu {
    queue = '<gpu-queue>'
    clusterOptions = '--gpus 1'
  }
}
```

:::tip
More information on how to configure Slurm in Nextflow can be found [here](https://www.nextflow.io/docs/latest/executor.html#slurm).
Depending on your cluster configuration, you might need to adjust the `clusterOptions` to one of the following:

- `--gpus 1` (as in the example above)
- `--gpus-per-node=1`
- `--gres=gpu:1`

:::

:::tip
If your jobs get assigned to the correct nodes, but the GPU is not utilized, you might need to add the following configuration:
`singularity.runOptions = '--no-mount tmp --writable-tmpfs --nv --env CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES --env ROCR_VISIBLE_DEVICES=$ROCR_VISIBLE_DEVICES --env ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK --env NVIDIA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES`

The first part (`--no-mount tmp --writable-tmpfs --nv`) is set by default in the `gpu` profile.
The rest of this configuration is needed in some cases to make the GPU visible to the container.
:::

For different executors, the configuration might look different.
Once a wider range of users have tested the GPU support, we will provide more detailed instructions for different executors.

### Ambient RNA correction

Ambient RNA correction removes contaminating RNA from cell-free droplets that can confound single-cell analysis.
The pipeline supports multiple ambient RNA correction methods that can be configured both globally and per-sample.

The pipeline allows you to select an ambient RNA correction method globally using the `--ambient_correction` parameter.
Available methods include `soupx` (default), `decontx`, `cellbender`, `scar`, or `none` to skip correction entirely.
SoupX requires an unfiltered matrix for each sample where ambient correction is enabled. For filtered-only samples, disable correction in the samplesheet or use `--ambient_correction decontx`.

> [!WARNING]
> If nf-core/scrnaseq already ran CellBender and you also enable downstream ambient correction on the same count matrix, you may apply two correction steps. Inspect the upstream outputs and disable one stage when appropriate.

```bash
nextflow run nf-core/scdownstream --ambient_correction decontx --input samplesheet.csv --outdir results
```

For finer control, you can disable ambient RNA correction for specific samples by setting `ambient_correction` to `false` in your samplesheet:

```csv title="samplesheet.csv"
sample,filtered,unfiltered,ambient_correction
sample1,/path/to/sample1_filtered.h5ad,/path/to/sample1.h5ad,true
sample2,/path/to/sample2_filtered.h5ad,/path/to/sample2.h5ad,false
```

By default, the pipeline stores ambient-corrected counts as additional layers in the AnnData object (e.g., `ambient_corrected_soupx`) while keeping the original raw counts in the `X` layer.
This means all downstream analysis including integration uses the raw counts, with corrected counts available for optional inspection.

If you want to use the ambient-corrected counts for integration instead, you can enable this behavior globally or per sample:

```bash
nextflow run nf-core/scdownstream --ambient_corrected_integration true --input samplesheet.csv --outdir results
```

```csv title="samplesheet.csv"
sample,filtered,unfiltered,ambient_corrected_integration
sample1,/path/to/sample1_filtered.h5ad,/path/to/sample1.h5ad,true
sample2,/path/to/sample2_filtered.h5ad,/path/to/sample2.h5ad,false
```

:::warning
When `ambient_corrected_integration` is enabled, the corrected counts replace the raw counts in the `X` layer, and the original raw counts are no longer available.
:::

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version.
When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since.
To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scdownstream
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data.
This ensures that a specific version of the pipeline code and software are used when you run your pipeline.
If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scdownstream releases page](https://github.com/nf-core/scdownstream/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.
Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.
For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile.
Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`.
This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `python_only`
  - Swaps R-based QC defaults for Python tools: scAR (`--ambient_correction`) instead of decontX, Scrublet (`--doublet_detection`) instead of scDblFinder, and scanpy HVGs (`--feature_selection hvgs`) instead of deviance feature selection
  - Combine with a software profile, e.g. `-profile docker,python_only`
  - scAR requires filtered and unfiltered matrices; use `--ambient_correction none` for filtered-only samples
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline.
Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
For input to be considered the same, not only the names must be identical but the files' contents as well.
For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`.
Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).
See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests.
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time.
For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original).
If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool.
By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects.
However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline.
Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository.
Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter.
You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs.
The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session.
The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
