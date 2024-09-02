include { SCVITOOLS_SCARCHES } from '../../modules/local/scvitools/scarches'

workflow TRANSFER {
    take:
    ch_transfer
    ch_base
    ch_inner
    
    main:
    ch_versions      = Channel.empty()
    ch_integrations  = Channel.empty()
    ch_obsm          = Channel.empty()

    SCVITOOLS_SCARCHES(
        ch_transfer.map{ meta, h5ad -> [[id: "scarches"], h5ad] },
        ch_base.map{ meta, h5ad, scvi_model -> [meta, scvi_model] },
        ch_inner
    )
    ch_versions      = ch_versions.mix(SCVITOOLS_SCARCHES.out.versions)
    ch_obsm          = ch_obsm.mix(SCVITOOLS_SCARCHES.out.obsm)
    ch_integrations  = ch_integrations.mix(SCVITOOLS_SCARCHES.out.h5ad)

    emit:
    obsm             = ch_obsm
    integrations     = ch_integrations

    versions         = ch_versions
}