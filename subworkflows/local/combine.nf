include { ADATA_MERGE      } from '../../modules/local/adata/merge'
include { ADATA_UPSETGENES } from '../../modules/local/adata/upsetgenes'
include { INTEGRATE        } from './integrate'

workflow COMBINE {

    take:
    ch_h5ad // channel: [ val(meta), h5ad ]

    main:

    ch_versions = Channel.empty()
    ch_obs = Channel.empty()
    ch_obsm = Channel.empty()
    ch_layers = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ADATA_MERGE(ch_h5ad.map { meta, h5ad -> h5ad }.collect()
                .map{ h5ads -> [[id: "merged"], h5ads] })
    ch_inner = ADATA_MERGE.out.inner
    ch_outer = ADATA_MERGE.out.outer
    ch_versions = ch_versions.mix(ADATA_MERGE.out.versions)

    ADATA_UPSETGENES(ch_outer)
    ch_versions = ch_versions.mix(ADATA_UPSETGENES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ADATA_UPSETGENES.out.multiqc_files)

    INTEGRATE(ch_inner)
    ch_versions = ch_versions.mix(INTEGRATE.out.versions)
    ch_integrations = INTEGRATE.out.integrations
    ch_obs = ch_obs.mix(INTEGRATE.out.obs)
    ch_obsm = ch_obsm.mix(INTEGRATE.out.obsm)

    emit:
    h5ad          = ch_outer
    integrations  = ch_integrations
    obs           = ch_obs
    obsm          = ch_obsm
    layers        = ch_layers

    multiqc_files = ch_multiqc_files
    versions      = ch_versions        // channel: [ versions.yml ]
}
