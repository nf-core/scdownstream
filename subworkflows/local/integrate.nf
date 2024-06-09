include { SCVITOOLS_SCVI } from '../../modules/local/scvitools/scvi'

workflow INTEGRATE {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()
    ch_integrations = Channel.empty()

    methods = params.integration_methods.split(',').collect{it.trim().toLowerCase()}

    if (methods.contains('scvi') || methods.contains('scanvi')) {
        SCVITOOLS_SCVI(ch_h5ad)
        ch_versions = ch_versions.mix(SCVITOOLS_SCVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCVI.out.h5ad
            .map{meta, h5ad -> [meta + [integration: 'scvi'], h5ad]})
    }

    emit:
    integrations = ch_integrations

    versions = ch_versions
}
