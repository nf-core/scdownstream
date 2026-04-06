include { UNTAR                } from '../../../modules/nf-core/untar'
include { SCIMILARITY_EMBED    } from '../../../modules/local/scimilarity/embed'
include { SCIMILARITY_ANNOTATE } from '../../../modules/local/scimilarity/annotate'

workflow SCIMILARITY {
    take:
    ch_h5ad           // channel: [ merged, h5ad ]
    scimilarity_model // channel: [ model ]

    main:
    ch_versions = channel.empty()
    ch_integrations = channel.empty()
    ch_obsm = channel.empty()
    ch_obs = channel.empty()

    ch_scimilarity_model = channel.value(
        [
            [id: 'scimilarity_model'],
            file(scimilarity_model),
        ]
    )

    if (!scimilarity_model) {
        error("scimilarity_model is required for scimilarity integration")
    }

    if (scimilarity_model.endsWith('.tar.gz')) {
        UNTAR (
            ch_scimilarity_model
        )
        ch_scimilarity_model = UNTAR.out.untar
    }

    SCIMILARITY_EMBED (
        ch_h5ad,
        ch_scimilarity_model,
    )
    ch_versions = ch_versions.mix(SCIMILARITY_EMBED.out.versions)
    ch_integrations = ch_integrations.mix(SCIMILARITY_EMBED.out.h5ad)
    ch_obsm = ch_obsm.mix(SCIMILARITY_EMBED.out.obsm)

    SCIMILARITY_ANNOTATE (
        SCIMILARITY_EMBED.out.h5ad,
        ch_scimilarity_model,
    )
    ch_versions = ch_versions.mix(SCIMILARITY_ANNOTATE.out.versions)
    ch_obs = ch_obs.mix(SCIMILARITY_ANNOTATE.out.obs)

    emit:
    integrations = ch_integrations // channel: [ integration, h5ad ]
    obs          = ch_obs          // channel: [ pkl ]
    obsm         = ch_obsm         // channel: [ pkl ]
    versions     = ch_versions     // channel: [ versions.yml ]
}
