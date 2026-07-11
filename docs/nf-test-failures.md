# nf-test failures overview

Summary of failing local module and subworkflow tests from parallel `nftu` runs on 2026-07-11.

**Profile:** `+apptainer,daisybio,keep_work`

## Log directories

| Run | Directory |
| --- | --------- |
| Initial module sweep (60 tests) | `.nf-test-parallel-20260711104815/` |
| Initial subworkflow sweep (20 tests) | `.nf-test-subworkflows-20260711121940/` |
| Rerun of previously failing tests (11 tests) | `.nf-test-rerun-20260711134343/` |

## Rerun summary (2026-07-11, 13:43)

Re-ran all 11 test files that failed in the initial sweeps.

| Result | Count |
| ------ | ----- |
| Fully passing | **4** |
| Still failing | **7** |

### Resolved since initial run

| Component | Initial failure | Rerun result |
| --------- | --------------- | ------------ |
| `modules/local/celltypist` | Stub test failed | All 4 tests pass |
| `modules/local/custom/doubletremoval` | `n_obs` assertion (6705 vs 7373) | All 3 tests pass |
| `subworkflows/local/doublet_detection` | `n_obs` assertion (6705 vs 7373) | All 4 tests pass |
| `subworkflows/local/sub_integrate` | Channel count mismatch (16 vs 17) | All 5 tests pass; 2 snapshots updated |

### Partially improved

| Component | Initial | Rerun |
| --------- | ------- | ----- |
| `subworkflows/local/combine` | 0/2 pass (channel mismatch) | 1/2 pass; stub fixed, full run still fails |
| `subworkflows/local/rank_genes_groups` | 0/3 pass | 1/3 pass; leiden test fixed |
| `subworkflows/local/per_group` | 1/13 fail (`full run`) | 1/13 fail (`full run - stub`); `full run` now passes |

### Still failing (7 test files)

| Component | Tests run | Failed | Snapshot updated |
| --------- | --------- | ------ | ---------------- |
| `modules/local/scanpy/cellcycle` | 2 | 1 | No |
| `modules/local/scran/normalization` | 2 | 1 | No |
| `subworkflows/local/combine` | 2 | 1 | Yes (partial) |
| `subworkflows/local/rank_genes_groups` | 3 | 2 | No |
| `subworkflows/local/quality_control` | 8 | 2 | No |
| `subworkflows/local/per_group` | 13 | 1 | No |
| `subworkflows/local/integrate` | 15 | 1 | No |

---

## Current failures (after rerun)

### Local modules

#### `scanpy/cellcycle`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run without failures - human | Failed | `SCANPY_CELLCYCLE` exit 1 |

Root cause (unchanged):

```
ValueError: zip() argument 2 is longer than argument 1
```

at `sc.pl.violin(adata, ["S_score", "G2M_score"], groupby="phase", ...)`. Cell-cycle genes are not found in `var_names`, so scoring/plotting fails.

Stub test passes.

#### `scran/normalization`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run without failures | Failed | `SCRAN_NORMALIZATION` exit 1 |

Root cause changed since initial run. Container now loads `scran`, but fails with:

```
Error in .local(x, ...) : size factors should be positive
```

at `logNormCounts` / `normalizeCounts`. Stub test passes.

### Local subworkflows

#### `combine`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run without failures - without base - stub | Passed | — |
| Should run without failures - without base | Failed | `COMBINE:INTEGRATE:FEATURE_SELECTION:SCRY_DEVIANCE (merged)` exit 1 |

Root cause:

```
Error: Selected deviant genes missing from input AnnData: CD74, IGHM, H1-4, ...
```

Gene symbols selected by deviance HVG selection are absent from the merged AnnData `var_names`.

#### `rank_genes_groups`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run differential expression with leiden clustering | Passed | — |
| Should run differential expression with label column | Failed | `assert workflow.out.uns.size() >= 1` — all output channels empty |
| Should run multiple rank-genes-groups methods in parallel | Failed | `Failed to load children for group '//' at address '96'` |

The leiden test that previously failed with empty dotplot data now passes. Remaining failures are an empty `uns` output for the label-column test and an HDF5 read error in the parallel-methods test.

#### `quality_control`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should apply sc-best-practices QC filtering when samplesheet defaults are present in meta | Failed | `assert adata.n_obs < 6000` — got **12663** cells |
| Should run with cell cycle scoring | Failed | `QUALITY_CONTROL:SCANPY_CELLCYCLE` exit 1 |

Cell-cycle failure is the same underlying issue as `scanpy/cellcycle` (violin plot / missing cell-cycle genes). Six other tests in this file pass.

#### `per_group`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run without failures - full run - stub | Failed | `Failed to open file '.../with_label_and_condition_pca.h5ad'. Is it a HDF5 file?` |

Previously failing `full run` test now passes. New failure is in the stub variant, likely an `nft-anndata` cache / stub fixture issue.

#### `integrate`

| Test | Result | Error |
| ---- | ------ | ----- |
| Should run without failures - deviance | Failed | `INTEGRATE:FEATURE_SELECTION:SCRY_DEVIANCE (test)` exit 1 |

Root cause changed since initial run (was `SCANPY_FILTER` / empty matrix). Now:

```
Error: Selected deviant genes missing from input AnnData: CD74, IGHM, H1-4, ...
```

Same deviance gene-matching issue as `combine`. Fourteen other tests pass, including `deviance - stub`.

---

## Initial sweep summary (for reference)

First parallel `nftu` runs before the rerun.

| Scope | Total | Passed | Failed |
| ----- | ----- | ------ | ------ |
| Local modules | 60 | 56 | 4 |
| Local subworkflows | 20 | 13 | 7 |

Of the 11 initially failing test files, 4 are now fully green and 3 are partially improved (see rerun summary above).

---

## Failure themes

| Theme | Affected components | Status | Suggested fix direction |
| ----- | ------------------- | ------ | ----------------------- |
| Cell-cycle scoring/plotting | `scanpy/cellcycle`, `quality_control` | Open | Fix gene symbol matching; handle empty phase groups in violin plot |
| Deviance HVG gene mismatch | `combine`, `integrate` | Open | Align deviance-selected gene symbols with AnnData `var_names` after merge/normalisation |
| SCRAN size factors | `scran/normalization` | Open | Investigate zero/negative size factors from `computeSumFactors` on test input |
| Rank genes empty output | `rank_genes_groups` (label column) | Open | Debug why `uns` channel is empty for label-based DE |
| HDF5 / cache read errors | `rank_genes_groups` (parallel), `per_group` (stub) | Open | Check stub fixtures and `nft-anndata` cache paths |
| QC cell-count assertion | `quality_control` | Open | Update expected `n_obs` threshold or fix sc-best-practices filtering defaults |
| Workflow input channel mismatch | `combine`, `sub_integrate` | **Resolved** | Fixed between initial run and rerun |
| Cell count assertion drift | `custom/doubletremoval`, `doublet_detection` | **Resolved** | Passing on rerun |
| Celltypist stub | `celltypist` | **Resolved** | Passing on rerun |
| Rank genes dotplot on empty results | `rank_genes_groups` (leiden) | **Resolved** | Leiden test passes on rerun |

---

## Re-running failed tests

```bash
cd /nfs/data/COST_IBD/scdownstream-dev
source ~/.bashrc

# All currently failing test files
nftu modules/local/scanpy/cellcycle/tests/main.nf.test \
      modules/local/scran/normalization/tests/main.nf.test \
      subworkflows/local/combine/tests/main.nf.test \
      subworkflows/local/rank_genes_groups/tests/main.nf.test \
      subworkflows/local/quality_control/tests/main.nf.test \
      subworkflows/local/per_group/tests/main.nf.test \
      subworkflows/local/integrate/tests/main.nf.test

# Single test file (no snapshot update)
nft <path/to/tests/main.nf.test>
```
