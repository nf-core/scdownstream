include { SCVITOOLS_REFERENCEMAPPING } from '../../modules/local/scvitools/referencemapping'

workflow TRANSFER {
    take:
    ch_transfer
    ch_base
    ch_inner

    main:
    ch_versions      = Channel.empty()
    ch_integrations  = Channel.empty()
    ch_obsm          = Channel.empty()

    SCVITOOLS_REFERENCEMAPPING(
        ch_transfer.map{ meta, h5ad -> h5ad }
            .combine(ch_base.map{ meta, h5ad, scvi_model, model_type -> [scvi_model, model_type] })
            .map{ h5ad, model, model_type -> [[id: model_type], h5ad, model, model_type]},
        ch_inner
    )
    ch_versions      = ch_versions.mix(SCVITOOLS_REFERENCEMAPPING.out.versions)
    ch_obsm          = ch_obsm.mix(SCVITOOLS_REFERENCEMAPPING.out.obsm)
    ch_integrations  = ch_integrations.mix(SCVITOOLS_REFERENCEMAPPING.out.h5ad)

    emit:
    obsm             = ch_obsm
    integrations     = ch_integrations

    versions         = ch_versions
}
