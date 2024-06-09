include { SCVITOOLS_SCVI   } from '../../modules/local/scvitools/scvi'
include { SCVITOOLS_SCANVI } from '../../modules/local/scvitools/scanvi'

workflow INTEGRATE {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()
    ch_obs = Channel.empty()
    ch_integrations = Channel.empty()

    methods = params.integration_methods.split(',').collect{it.trim().toLowerCase()}

    if (methods.contains('scvi') || methods.contains('scanvi')) {
        SCVITOOLS_SCVI(ch_h5ad)
        ch_versions = ch_versions.mix(SCVITOOLS_SCVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCVI.out.h5ad
            .map{meta, h5ad -> [[id: 'scvi'], h5ad]})

        if (methods.contains('scanvi')) {
            SCVITOOLS_SCANVI(ch_h5ad, SCVITOOLS_SCVI.out.model.collect())
            ch_versions = ch_versions.mix(SCVITOOLS_SCANVI.out.versions)
            ch_integrations = ch_integrations.mix(SCVITOOLS_SCANVI.out.h5ad
                .map{meta, h5ad -> [[id: 'scanvi'], h5ad]})
            ch_obs = ch_obs.mix(SCVITOOLS_SCANVI.out.obs)
        }
    }

    emit:
    integrations = ch_integrations
    obs = ch_obs

    versions = ch_versions
}
