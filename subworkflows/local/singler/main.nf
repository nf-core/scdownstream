include { CELLDEX_FETCHREFERENCE } from '../../../modules/local/celldex/fetchreference'
include { CELLTYPES_SINGLER      } from '../../../modules/local/singler'

workflow SINGLER {
    take:
    ch_h5ad      // channel: [ meta, h5ad, symbol_col, counts_layer ]
    ch_reference // channel: [ meta, reference ]

    main:
    ch_obs = channel.empty()

    ch_reference = ch_reference.branch { _meta, ref ->
        files: file(ref).exists() && file(ref).isFile()
        names: true
    }

    CELLDEX_FETCHREFERENCE (
        ch_reference.names
            .map { meta, ref ->
                if (!meta.version) {
                    error "If you specify a celldex reference, you also need to specify a version"
                }
                [meta, ref, meta.version]
            }
    )

    // Bring the branches back together
    ch_reference = ch_reference.files.mix(CELLDEX_FETCHREFERENCE.out.tar)

    CELLTYPES_SINGLER (
        ch_h5ad,
        ch_reference
            .map {
                meta, ref ->
                [[id: "singler"], meta.id, meta.label, ref]
            }
            .groupTuple()
            .collect()
    )
    ch_obs = ch_obs.mix(CELLTYPES_SINGLER.out.obs)

    emit:
    obs      = ch_obs
}
