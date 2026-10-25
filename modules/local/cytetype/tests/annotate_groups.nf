include { ADATA_SPLITCOL                        } from '../../adata/splitcol/main'
include { ADATA_EXTEND as ADD_GROUPS            } from '../../adata/extend/main'
include { anndata                               } from 'plugin/nf-anndata'
include { SCANPY_RANKGENESGROUPS as RANK_GLOBAL } from '../../scanpy/rankgenesgroups/main'
include { SCANPY_RANKGENESGROUPS as RANK_SUBSET } from '../../scanpy/rankgenesgroups/main'
include { CYTETYPE as CYTETYPE_GLOBAL           } from '../main'
include { CYTETYPE as CYTETYPE_SUBSET           } from '../main'

workflow ANNOTATE_GROUPS {
    take:
    h5ad

    main:
    group_obs = h5ad
        .map { _meta, input ->
            def ad = anndata(input)
            def labels = ad.obs.get('leiden').data
            def names = ad.obs.index.values
            def rows = (0..<ad.n_obs).collect { index ->
                "${names[index]},${labels[index]},${labels[index]}"
            }
            'cell,pca-global-0.5_leiden,pca-subset_with_underscore-0.5_leiden\n' + rows.join('\n') + '\n'
        }
        .collectFile(name: 'groups.csv')
    ADD_GROUPS(
        h5ad.combine(group_obs).map { meta, input, obs ->
            [meta + [id: 'grouped'], input, [obs], [], [], [], [], []]
        }
    )
    ADATA_SPLITCOL(ADD_GROUPS.out.h5ad, 'sample')
    subset = ADATA_SPLITCOL.out.h5ad.map { meta, files ->
        [meta + [id: 'subset'], files.sort { item -> item.name }.first()]
    }
    RANK_GLOBAL(ADD_GROUPS.out.h5ad, 'pca-global-0.5_leiden', [[], []], 'wilcoxon', 'rank_global')
    RANK_SUBSET(subset, 'pca-subset_with_underscore-0.5_leiden', [[], []], 'wilcoxon', 'rank_subset')
    CYTETYPE_GLOBAL(RANK_GLOBAL.out.h5ad, 'index', 'Offline global annotation', 'pca-global-0.5_leiden', 'rank_global', 'pca', '0.5')
    CYTETYPE_SUBSET(RANK_SUBSET.out.h5ad, 'index', 'Offline subset annotation', 'pca-subset_with_underscore-0.5_leiden', 'rank_subset', 'pca', '0.5')

    emit:
    obs         = CYTETYPE_GLOBAL.out.obs.mix(CYTETYPE_SUBSET.out.obs)
    subset_h5ad = CYTETYPE_SUBSET.out.h5ad
    base_index  = h5ad.map { _meta, input -> [anndata(input).obs.index.values.toList()] }
}
