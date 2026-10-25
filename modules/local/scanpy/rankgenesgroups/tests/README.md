# Expression-source regression

`explicit_x.nf.test` compares real process outputs from the D01 matched fixture pair. Both files have identical counts in `X`, ordered axes and grouping metadata. The second file has a `.raw` marker pattern that reverses the exclusive markers in `X`.

The fixtures and their [generator recipe](https://raw.githubusercontent.com/nictru/test-datasets/55565654bd8b77673cc6d8ebab39a16aea1e6822/anndata-variations/DE_RAW_FIXTURES.md) are pinned to fork commit `55565654bd8b77673cc6d8ebab39a16aea1e6822`. After the dataset PR merges, replace these URLs with the immutable upstream commit containing the same files, without changing the global test-data base.

Wilcoxon and t-test runs must produce identical standard DE rows, filtered marker CSVs and output schemas for both inputs. Both rank and dot plots must be emitted. Equality is exact because the selected input matrices and method parameters are identical. The existing logreg test also checks that unavailable fold changes and pvalues remain null in the standard table.

The pinned Scanpy 1.12.3 `rank_genes_groups` and `filter_rank_genes_groups` APIs accept `use_raw=False`. The dotplot wrapper forwards this option to the plotting API. Explicit selection preserves `.raw` and retained layers, while ranking, marker filtering and dotplot expression all use transformed `X`.

The current nft-anndata metadata interface exposes `uns` keys but not nested rank parameters. The automated test checks the rank key and ordered observation identifiers; it cannot inspect `uns['rank_genes_groups']['params']['use_raw']`. Scratch inspection of all four real matched-pair outputs confirmed `use_raw=False` for both methods. The conflicting-input outputs retained `.raw`; the absent-input outputs did not acquire it. No HDF5 reader or test-only assertion process is added for this automation gap.
