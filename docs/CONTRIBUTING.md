---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/scdownstream`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/scdownstream.

This page describes the recommended nf-core way to contribute to both nf-core/scdownstream and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/scdownstream.

> [!NOTE]
> If you need help using or modifying nf-core/scdownstream, ask on the nf-core Slack [#scdownstream](https://nfcore.slack.com/channels/scdownstream) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/scdownstream/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/scdownstream repository](https://github.com/nf-core/scdownstream) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focused as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/scdownstream documentation](https://nf-co.re/scdownstream/usage) or ask on the nf-core [#scdownstream](https://nfcore.slack.com/channels/scdownstream) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/scdownstream without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/scdownstream/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/scdownstream code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

### Harshil alignment

Use [Harshil alignment](https://nf-co.re/docs/developing/documentation/harshil-alignment) where applicable.
Do not add more spaces than the widest element requires.

### Linting and formatting

Stage your changes and run [`prek`](https://github.com/j178/prek) before each commit.
The hooks in `.pre-commit-config.yaml` run Prettier, `nextflow lint`, and `ruff` for Python code.
The `ruff` rules are configured in `ruff.toml`.

### Local modules

Local modules follow the [nf-core module specifications](https://nf-co.re/docs/specifications/components/modules/general).

### Meta map fields

In addition to `id`, local modules and subworkflows use the following meta fields:

| Field                             | Set by                                               | Meaning                                                                           |
| --------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| `batch_col`                       | Samplesheet                                          | `obs` column holding batch labels, used for batch-aware QC and feature selection. |
| `symbol_col`                      | Samplesheet                                          | `var` column holding gene symbols.                                                |
| `counts_layer`                    | Samplesheet                                          | Layer holding raw counts. Defaults to `X`.                                        |
| `condition_col`, `donor_col`      | Samplesheet                                          | `obs` columns used for differential expression and pseudobulking.                 |
| `integration`                     | Integration step                                     | Name of the integration method that produced the embedding.                       |
| `subset`                          | `CLUSTER_TARGETS`, `SUB_INTEGRATE`                   | `global` for the full dataset, otherwise the name of the subset.                  |
| `resolution`                      | Clustering step                                      | Leiden resolution of the clustering.                                              |
| `obs_key`                         | `workflows/scdownstream.nf`                          | `obs` column that defines the groups for per-group analyses.                      |
| `analyses`, `de_methods_resolved` | Analysis plan (`utils_nfcore_scdownstream_pipeline`) | Per-group analyses and DE methods to run. `null` means all.                       |

Channel operations must preserve these fields.
Add new fields to this table when you introduce them.

### Process templates

Many processes in this pipeline use templates containing Python or R scripts.
Nextflow interprets strings prefixed with `$` as Nextflow variables - this is especially relevant for R scripts.
The process fails if such a variable is not defined in the associated Nextflow process.
To use a literal `$` character in a template, escape it with a backslash (`\$`).
In R templates, this applies to every use of the `$` operator, for example `adata\$as_SingleCellExperiment()` or `rowData(sce)\$highly_deviant`.

Environment variables that apply to all Python processes, such as `NUMBA_CACHE_DIR` and `MPLCONFIGDIR`, are set in the `env` scope of `nextflow.config`.
Do not set them in individual templates.
Pass `task.cpus` to the tool explicitly, for example with `threadpool_limits(int("${task.cpus}"))` in Python.
Set a fixed random seed for stochastic methods, for example with `set.seed()` in R.

Templates should follow this structure:

1. Read inputs.
2. Perform processing.
3. Write outputs:
   - Primary outputs, such as updated AnnData files.
   - Fragments, such as added `obs` columns or layers, in relevant formats.
   - MultiQC inputs.
   - `versions.yml` containing the versions of the language and all relevant packages.

### Fragments

A fragment contains only the data that a module adds to an AnnData object.
`workflows/scdownstream.nf` collects fragments on the `ch_obs`, `ch_var`, `ch_obsm`, `ch_obsp`, and `ch_uns` channels, and on their `_per_sample` counterparts (including `ch_layers_per_sample`) before combination.
`ADATA_EXTEND` merges them into the H5AD files, once per sample after quality control (`FINALIZE_QC_ANNDATAS`) and once at the end of the pipeline (`FINALIZE`).
Emit each fragment on a channel named after its AnnData slot and use a format that `modules/local/adata/extend` can read:

| Slot     | Format                                               | Index                                 |
| -------- | ---------------------------------------------------- | ------------------------------------- |
| `obs`    | Pickled or CSV `pandas.DataFrame`                    | `obs_names`                           |
| `var`    | Pickled or CSV `pandas.DataFrame`                    | `var_names`                           |
| `obsm`   | Pickled `pandas.DataFrame`, file name `X_<name>.pkl` | `obs_names`                           |
| `obsp`   | NumPy `.npy` file containing a sparse matrix         | None                                  |
| `uns`    | Pickled Python object                                | None                                  |
| `layers` | `.npz`, `.mtx`, or `.npy` matrix                     | Same cell and gene order as the input |

The file stem becomes the key in `obsm`, `obsp`, `uns`, and `layers`.

### MultiQC custom content

Modules that produce plots for the MultiQC report write a `<prefix>_mqc.json` file with an embedded base64 PNG image.
Set `parent_id`, `parent_name`, and `parent_description` from `meta.integration` so that MultiQC groups the section with the other results of the same integration.
`modules/local/scanpy/leiden/templates/leiden.py` contains a reference implementation.
Add new sections to `assets/multiqc_config.yml` in the order of the pipeline steps.

### Output directories

`conf/modules.config` publishes results into numbered top-level directories that follow the order of the pipeline steps, for example `02_quality_control` or `05_cluster_dimred`.
Place outputs of new steps into the matching directory.
Guard intermediate outputs with `enabled: params.save_intermediates`.
Update `docs/output.md` with the new files.

### Python-only profile

The `python_only` profile in `conf/python_only.config` replaces R-based default tools with Python alternatives.
If you add an R-based tool as a new default, add a Python alternative to this profile.
If no alternative exists, document the limitation in the comment block of `conf/python_only.config`.

### Testing conventions

CI runs two nf-test tiers:

- `modules_local,subworkflows_local` for local module and subworkflow tests.
- `pipeline` for end-to-end tests in `tests/`.

Tag new tests accordingly.

[`docs/reproducibility.md`](reproducibility.md) classifies every local module and subworkflow by the reproducibility of its outputs and assigns a snapshot strategy to each.
When you add a module or subworkflow, add a row to this file and write the test according to the assigned strategy.
Snapshot files must never be edited by hand; regenerate them with `nf-test test --update-snapshot`.
