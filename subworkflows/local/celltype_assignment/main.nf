include { samplesheetToList    } from 'plugin/nf-schema'
include { SINGLER              } from '../singler'
include { CELLTYPES_CELLTYPIST } from '../../../modules/local/celltypes/celltypist'

workflow CELLTYPE_ASSIGNMENT {
    take:
    ch_h5ad           // channel: [ meta, h5ad, symbol_col, counts_layer ]
    celldex_reference //   value: string
    celltypist_model  //   value: string

    main:
    ch_versions = channel.empty()
    ch_obs = channel.empty()

    if (celldex_reference ) {
        SINGLER (
            ch_h5ad,
            channel.fromList(samplesheetToList(
                celldex_reference,
                "${projectDir}/assets/schema_singler.json")
            )
        )
        ch_obs = ch_obs.mix(SINGLER.out.obs)
        ch_versions = ch_versions.mix(SINGLER.out.versions)
    }

    if (celltypist_model) {
        celltypist_models = channel.value(celltypist_model
            .split(',')
            .collect{ it -> it.trim() }
        )

        CELLTYPES_CELLTYPIST (
            ch_h5ad.map { meta, h5ad, symbol_col, _counts_layer -> [meta, h5ad, symbol_col] },
            celltypist_models
        )
        ch_obs = ch_obs.mix(CELLTYPES_CELLTYPIST.out.obs)
        ch_versions = ch_versions.mix(CELLTYPES_CELLTYPIST.out.versions)
    }

    emit:
    obs      = ch_obs      // channel: [ meta, pkl ]
    versions = ch_versions // channel: [ versions.yml ]
}
