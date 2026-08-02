include { SCANPY_RANKGENESGROUPS   } from '../../../modules/local/scanpy/rankgenesgroups'
include { anndata                   } from 'plugin/nf-anndata'
include { rankGenesGroupsMethods    } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsMethodSlug } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsMethodKey  } from '../utils_nfcore_scdownstream_pipeline'
include { hasMultipleObsGroups      } from '../utils_nfcore_scdownstream_pipeline'

workflow RANK_GENES_GROUPS {
    take:
    ch_h5ad           // channel: [ meta, h5ad ], anndata objects with obs_key and condition_col in meta
    interesting_genes //   value: string (path) or []

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_settings = ch_h5ad.map { meta, h5ad ->
        def obs_key = meta.obs_key
        def condition_col = meta.condition_col
        def ad = anndata(h5ad)
        def rgg_methods = meta.de_methods_resolved.intersect(rankGenesGroupsMethods())

        def conditions = ad.obs[condition_col].unique().toList()
        def labels = ad.obs[obs_key].unique().toList()

        return [
            meta,
            h5ad,
            condition_col,
            conditions.size() > 1 ? conditions : [],
            obs_key,
            labels.size() > 1 ? labels : [],
            rgg_methods,
        ]
    }

    // Structure: [meta, h5ad, filter_col, filter_val, obs_key]
    ch_global_comparisons = ch_settings
        .map { meta, h5ad, _condition_col, _conditions, obs_key, _labels, _rgg_methods ->
            [meta, h5ad, '', '', obs_key]
        }

    ch_condition_labels = ch_settings.transpose(by: 3)
        .map { meta, h5ad, condition_col, condition, obs_key, _labels, _rgg_methods ->
            [meta, h5ad, condition_col, condition, obs_key]
        }

    ch_label_conditions = ch_settings.transpose(by: 5)
        .map { meta, h5ad, condition_col, _conditions, obs_key, label, _rgg_methods ->
            [meta, h5ad, obs_key, label, condition_col]
        }

    ch_filtered_comparisons = ch_label_conditions.mix(ch_condition_labels).map { meta, h5ad, filter_col, filter_val, obs_key ->
        [meta, h5ad, filter_col as String, filter_val as String, obs_key]
    }

    ch_comparisons = ch_global_comparisons.mix(ch_filtered_comparisons)

    ch_rankgenesgroups = ch_comparisons.filter { meta, h5ad, filter_col, filter_val, obs_key ->
            hasMultipleObsGroups(meta, h5ad, 'SCANPY_RANKGENESGROUPS', obs_key, filter_col, filter_val, 2)
        }
        .flatMap { meta, h5ad, filter_col, filter_val, obs_key ->
            def rgg_methods = meta.de_methods_resolved.intersect(rankGenesGroupsMethods())
            def single_method = rgg_methods.size() == 1
            def base_id = filter_col && filter_val ? "${obs_key}:${filter_col}:${filter_val}" : obs_key
            def comparison_scope = filter_col && filter_val ? 'filtered' : 'global'

            rgg_methods.collect { method ->
                def method_slug = rankGenesGroupsMethodSlug(method)
                def rank_key = rankGenesGroupsMethodKey(method, single_method)
                [
                    meta: meta + [
                        id: "${base_id}:${method_slug}",
                        comparison_scope: comparison_scope,
                        de_method: method,
                        rank_key: rank_key,
                    ],
                    h5ad: h5ad,
                    filter_col: filter_col,
                    filter_val: filter_val,
                    obs_key: obs_key,
                    method: method,
                    rank_key: rank_key,
                ]
            }
        }
        .multiMap { task ->
            h5ad: [task.meta, task.h5ad]
            obs_key: task.obs_key
            filter: [task.filter_col, task.filter_val]
            method: task.method
            rank_key: task.rank_key
        }

    SCANPY_RANKGENESGROUPS(
        ch_rankgenesgroups.h5ad,
        ch_rankgenesgroups.obs_key,
        ch_rankgenesgroups.filter,
        ch_rankgenesgroups.method,
        ch_rankgenesgroups.rank_key,
        interesting_genes ?: [],
    )
    ch_uns           = ch_uns.mix(SCANPY_RANKGENESGROUPS.out.uns.flatten())
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_RANKGENESGROUPS.out.multiqc_files.flatten())

    emit:
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    h5ad          = SCANPY_RANKGENESGROUPS.out.h5ad
                        .filter { meta, _h5ad -> meta.comparison_scope == 'global' }
                                     // channel: [ meta, h5ad ] — global comparisons only
}
