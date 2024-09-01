# nf-core/scdownstream: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/scdownstream/usage](https://nf-co.re/scdownstream/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Filtered and unfiltered matrices

Throughout this documentation, you will find references to `filtered` and `unfiltered` matrices.
The `unfiltered` matrices are matrices which still contain empty droplets, whereas the `filtered` matrices have been filtered for empty droplets. A more technical definition can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices). `CellRanger` provides you with both matrices, whereas other quantification tools only provide you with the `unfiltered` matrix.
The pipeline can handle the following cases:

1. You have both `filtered` and `unfiltered` matrices: Provide both matrices in the samplesheet and the pipeline will use the `unfiltered` matrix for ambient RNA removal and the `filtered` matrix for all other steps.
2. You only have the `filtered` matrix: Provide the `filtered` matrix in the samplesheet and the pipeline will use it for all steps. In this case, only `decontX` can be used for ambient RNA removal, as all other methods require the `unfiltered` matrix.
3. You only have the `unfiltered` matrix: Provide the `unfiltered` matrix in the samplesheet and the pipeline will automatically create a `filtered` matrix by identifying empty droplets using `CellBender`.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with at least 2 columns, and a header row as shown in the examples below.

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
sample,filtered,unfiltered,batch_col,label_col,unknown_label,min_genes,min_cells,min_counts_cell,min_counts_gene
sample1,/absolute/path/to/sample1_filtered.h5ad,/absolute/path/to/sample1.h5ad,batch,cell_type,unknown,1,2,3,4
sample2,relative/path/to/sample2_filtered.rds,relative/path/to/sample2.rds,batch_id,annotation,unannotated,5,6,7,8
sample3,/absolute/path/to/sample3_filtered.csv,/absolute/path/to/sample3.csv,,,,9,10,11,12
```

For CSV input files, specifying the `batch_col`, `label_col`, and `unknown_label` columns will not have any effect, as no additional metadata is available in the CSV file.

| Column            | Description                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`          | Unique sample identifier. Will be added to the pipeline output objects as `sample` column.                                                                                                                                                                                                                                                                                                                           |
| `filtered`        | May contain paths to `h5ad`, `h5`, `rds`, or `csv` files. `rds` files may contain any object that can be converted to a `SingleCellExperiment` using the [Seurat `as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.singlecellexperiment) function. `csv` files should contain a matrix with genes as columns and cells as rows.                                                                   |
| `unfiltered`      | Same as `file`, but for the unfiltered cellranger or nf-core/scrnaseq output. If not provided, only `decontX` can be used for ambient RNA removal.                                                                                                                                                                                                                                                                   |
| `batch_col`       | Column in the input file containing batch information. Defaults to `batch`. If the column does not exist in the input object, the pipeline will create a new column and put the sample identifier in it. If the `batch_col` is something else than `batch`, it will be renamed to `batch` during pipeline execution.                                                                                                 |
| `symbol_col`      | Column in the input file containing gene symbol information. Defaults to `index`. There are two special values that can be used: `index` and `none`. `index` will use the row names of the matrix as gene symbols. `none` will trigger the pipeline to perform gene symbol conversion (this is not supported yet). The values from `symbol_col` will be copied to a column `gene_symbols` during pipeline execution. |
| `label_col`       | Column in the input file containing cell type information. Defaults to `label`. If the column does not exist in the input object, the pipeline will create a new column and put `unknown` in it. If the `label_col` is something else than `label`, it will be renamed to `label` during pipeline execution.                                                                                                         |
| `unknown_label`   | Value in the `label_col` column that should be considered as unknown. Defaults to `unknown`. If the `unknown_label` is something else than `unknown`, it will be renamed to `unknown` during pipeline execution. If trying to perform integration with scANVI, more than one unique label other than `unknown` must exist in the input data.                                                                         |
| `min_genes`       | Minimum number of genes required for a cell to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                       |
| `min_cells`       | Minimum number of cells required for a gene to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                       |
| `min_counts_cell` | Minimum number of counts required for a cell to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                      |
| `min_counts_gene` | Minimum number of counts required for a gene to be considered. Defaults to `1`.                                                                                                                                                                                                                                                                                                                                      |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scdownstream --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scdownstream -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### GPU acceleration

:::warning{title="Experimental feature"}
This is an experimental feature and may produce errors. If you encounter any issues, please report them on the [nf-core/scdownstream GitHub repository](https://github.com/nf-core/scdownstream/issues/new?assignees=&labels=bug&projects=&template=bug_report.yml).
:::

:::info{title="Using conda"}
GPU acceleration is not available when using conda for dependency management.
:::

Tools with implemented support for GPU acceleration are:
- cellbender
- scVI
- scANVI
- scAR
- solo

To utilize GPU acceleration, you need to specify the `gpu` profile. This will make the tool steps use cuda-enabled environments and it will tell the tools to use the GPU. All processes which support GPU acceleration are marked with the `process_gpu` label.

You also need to make sure that the tasks are run on a machine with a GPU. If all tasks are run on a machine with a GPU, no further action is needed. If you are running the pipeline on a slurm cluster, where there is dedicated queue for GPU jobs, you need additional configuration that might look like this:

```bash
process {
  withLabel:process_gpu {
    queue = '<gpu-queue>'
    clusterOptions = '--gpus 1'
  }
}
```

For different executors, the configuration might look different. Once a wider range of users have tested the GPU support, we will provide more detailed instructions for different executors.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scdownstream
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scdownstream releases page](https://github.com/nf-core/scdownstream/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
