include { ADATA_MERGE           } from '../../modules/local/adata/merge'
include { ADATA_UPSETGENES      } from '../../modules/local/adata/upsetgenes'
include { INTEGRATE             } from './integrate'
include { ADATA_MERGEEMBEDDINGS } from '../../modules/local/adata/mergeembeddings'

workflow COMBINE {

    take:
    ch_h5ad // queue channel: [ val(meta), path(h5ad) ]
    ch_base  // value channel: [ val(meta), path(h5ad) ]
    ch_reference_model // value channel: [ val(meta), path(model) ]

    main:

    ch_versions      = Channel.empty()
    ch_obs           = Channel.empty()
    ch_obsm          = Channel.empty()
    ch_layers        = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ADATA_MERGE(
        ch_h5ad.map { meta, h5ad -> [[id: "merged"], h5ad] }.groupTuple(),
        ch_base
    )
    ch_outer         = ADATA_MERGE.out.outer
    ch_versions      = ch_versions.mix(ADATA_MERGE.out.versions)

    ADATA_UPSETGENES(ch_outer)
    ch_versions      = ch_versions.mix(ADATA_UPSETGENES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ADATA_UPSETGENES.out.multiqc_files)

    INTEGRATE( ADATA_MERGE.out.integrate, ch_reference_model )
    ch_versions      = ch_versions.mix(INTEGRATE.out.versions)

    if (params.base_adata) {
        ADATA_MERGEEMBEDDINGS(
            INTEGRATE.out.integrations
            .combine(
                ch_base.map{meta, base -> base}
            ).combine(
                ADATA_MERGE.out.inner.map{meta, inner -> inner}
            )
        )
        ch_versions      = ch_versions.mix(ADATA_MERGEEMBEDDINGS.out.versions)
        ch_integrations  = ADATA_MERGEEMBEDDINGS.out.h5ad
        ch_obs           = ch_obs.mix(ADATA_MERGEEMBEDDINGS.out.obs)
        ch_obsm          = ch_obsm.mix(ADATA_MERGEEMBEDDINGS.out.obsm)
    } else {
        ch_integrations  = INTEGRATE.out.integrations
        ch_obs           = ch_obs.mix(INTEGRATE.out.obs)
        ch_obsm          = ch_obsm.mix(INTEGRATE.out.obsm)
    }

    ch_integrations = ch_integrations
        .map{meta, file -> [meta + [integration: meta.id], file]}

    emit:
    h5ad             = ch_outer
    h5ad_inner       = ADATA_MERGE.out.inner
    integrations     = ch_integrations
    obs              = ch_obs
    obsm             = ch_obsm
    layers           = ch_layers

    multiqc_files    = ch_multiqc_files
    versions         = ch_versions        // channel: [ versions.yml ]
}
