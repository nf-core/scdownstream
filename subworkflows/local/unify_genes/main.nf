include { HUGOUNIFIER_GET   } from '../../../modules/local/hugounifier/get'
include { HUGOUNIFIER_APPLY } from '../../../modules/local/hugounifier/apply'

workflow UNIFY_GENES {
    take:
    ch_h5ad // channel: [ meta, h5ad ]

    main:

    HUGOUNIFIER_GET(
        ch_h5ad
            .map { meta, h5ad -> [[id: 'hugo-unifier'], meta.id, h5ad] }
            .groupTuple()
    )

    ch_changes = HUGOUNIFIER_GET.out.changes
        .map { _meta, changes -> [changes] }
        .flatten()
        // Extract ID from file name
        .map { changes -> [changes.baseName, changes] }
        // Join with AnnData based on ID
        .join(
            ch_h5ad.map{ meta, h5ad -> [meta.id, meta, h5ad] }
        )
        // Prettify output
        .map { _id, changes, meta, h5ad ->
            [meta, h5ad, changes]
        }

    HUGOUNIFIER_APPLY (
        ch_changes
    )
    ch_h5ad = HUGOUNIFIER_APPLY.out.h5ad

    emit:
    h5ad = ch_h5ad // channel: [ meta, h5ad ]
}
