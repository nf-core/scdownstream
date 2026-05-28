# nf-core/test-datasets update — `scdownstream/extension_base`

Copy the contents of `extension_base/` into the **`scdownstream` branch** of [nf-core/test-datasets](https://github.com/nf-core/test-datasets), replacing the existing files in `scdownstream/extension_base/`.

## Files

| File                                    | Purpose                                                      |
| --------------------------------------- | ------------------------------------------------------------ |
| `extension_base/model.pt`               | scVI checkpoint for reference mapping / extension            |
| `extension_base/merged.h5ad`            | Finalized atlas (`base_adata`) for extension                 |
| `extension_base/harmony_reference.h5ad` | Symphony reference for Harmony reference mapping / extension |

All three must come from the **same pipeline run** (see below).

> **Note:** If `extension_base/` already contains files but you have not run `collect-artifacts.sh` yet, `merged.h5ad` and `model.pt` may still be the current test-datasets versions and `harmony_reference.h5ad` from a Harmony-only build. **Do not open the test-datasets PR until you have run a unified build and `collect-artifacts.sh`** — all three files must be replaced together.

## How these were generated

```bash
# From the scdownstream repo root, with nf-core conda env active:
nextflow run main.nf -profile test,apptainer -params-file contrib/nf-core-test-datasets/build.params.json

# Populate extension_base/ from the build output:
./contrib/nf-core-test-datasets/collect-artifacts.sh
```

Build parameters match the consolidated pipeline tests (`main_pipeline_reference_mapping.nf.test`, `main_pipeline_extend.nf.test`):

- Input: `samplesheet.csv` (full atlas)
- Integration: `scvi,harmony`
- HVGs: 500

## PR checklist (test-datasets repo)

1. Check out branch `scdownstream`.
2. Ensure Git LFS is enabled (`git lfs install`).
3. Copy `extension_base/*` into `scdownstream/extension_base/` (overwrite `model.pt` and `merged.h5ad`, add `harmony_reference.h5ad`).
4. Add or extend `.gitattributes` on the test-datasets repo:

   ```
   scdownstream/extension_base/*.h5ad filter=lfs diff=lfs merge=lfs -text
   ```

5. Commit and open PR against `scdownstream`.
6. After merge, re-run `nftu` on the scdownstream pipeline reference-mapping and extend tests.

## Pipeline version

Record the scdownstream commit/tag used when generating these files in your PR description.
