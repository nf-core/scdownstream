include { ADATA_SPLITCOL as SPLITCOL } from '../../../modules/local/adata/splitcol'

workflow CLUSTER_TARGETS {
    take:
    ch_input    // channel: [ meta, h5ad ]
    per_label   //   value: boolean
    global      //   value: boolean
    split_col   //   value: string

    main:
    ch_targets = channel.empty()

    ch_input_by_subset = ch_input.branch { meta, _h5ad ->
        already_split: meta.subset != null
        needs_split: true
    }

    ch_targets = ch_targets.mix(ch_input_by_subset.already_split)

    if (global) {
        ch_targets = ch_targets.mix(
            ch_input_by_subset.needs_split
                .map { meta, h5ad -> [meta + [subset: 'global'], h5ad] }
        )
    }

    if (per_label) {
        SPLITCOL (
            ch_input_by_subset.needs_split,
            split_col
        )

        ch_targets = ch_targets.mix(
            SPLITCOL.out.h5ad
                .transpose()
                .map { meta, h5ad -> [meta + [subset: h5ad.simpleName], h5ad] }
        )
    }

    ch_targets = ch_targets.map { meta, h5ad ->
        def cluster_id = meta.subset != null
            ? meta.integration + '-' + meta.subset
            : meta.integration
        [meta + [id: cluster_id], h5ad]
    }

    emit:
    targets = ch_targets // channel: [ meta, h5ad ]
}
