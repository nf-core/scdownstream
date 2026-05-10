include { SCANPY_RANKGENESGROUPS } from '../../../modules/local/scanpy/rankgenesgroups'
include { anndata } from 'plugin/nf-anndata'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    ch_h5ad // channel: [ meta, h5ad ], anndata objects with obs_key and condition_col in meta

    main:
    ch_versions      = channel.empty()
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_settings = ch_h5ad.map { meta, h5ad ->
        def obs_key = meta.obs_key
        def condition_col = meta.condition_col
        def ad = anndata(h5ad)

        def conditions = ad.obs[condition_col].unique().toList()
        def labels = ad.obs[obs_key].unique().toList()

        return [
            meta,
            h5ad,
            condition_col,
            conditions.size() > 1 ? conditions : [],
            obs_key,
            labels.size() > 1 ? labels : []
        ]
    }

    // Structure: [meta, h5ad, filter_col, filter_val, obs_key]
    ch_comparisons = channel.empty()

    ch_global_labels = ch_settings
        .map { meta, h5ad, _condition_col, _conditions, obs_key, _labels ->
            [meta + [id: obs_key], h5ad, [], [], obs_key]
        }
    ch_comparisons = ch_comparisons.mix(ch_global_labels)

    ch_condition_labels = ch_settings.transpose(by: 3)
        .map { meta, h5ad, condition_col, condition, obs_key, _labels ->
            [meta, h5ad, condition_col, condition, obs_key]
        }

    ch_label_conditions = ch_settings.transpose(by: 5)
        .map { meta, h5ad, condition_col, _conditions, obs_key, label ->
            [meta, h5ad, obs_key, label, condition_col]
        }
    ch_comparisons = ch_comparisons.mix(
        ch_label_conditions.mix(ch_condition_labels).map { meta, h5ad, filter_col, filter_val, obs_key ->
            [meta + [id: obs_key + ":" + filter_col + ":" + filter_val], h5ad, filter_col, filter_val, obs_key]
        }
    )

    ch_rankgenesgroups = ch_comparisons.multiMap { meta, h5ad, filter_col, filter_val, obs_key ->
        h5ad: [meta, h5ad]
        obs_key: obs_key
        filter: [filter_col, filter_val]
    }

    SCANPY_RANKGENESGROUPS(
        ch_rankgenesgroups.h5ad,
        ch_rankgenesgroups.obs_key,
        ch_rankgenesgroups.filter,
        params.rankgenesgroups_method
    )
    ch_versions      = ch_versions.mix(SCANPY_RANKGENESGROUPS.out.versions)
    ch_uns           = ch_uns.mix(SCANPY_RANKGENESGROUPS.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_RANKGENESGROUPS.out.multiqc_files)

    emit:
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions      // channel: [ versions.yml ]
}
