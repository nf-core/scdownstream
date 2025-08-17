include { INTEGRATE             } from './integrate'
include { ADATA_MERGEEMBEDDINGS } from '../../modules/local/adata/mergeembeddings'
include { ADATA_MERGE           } from '../../modules/local/adata/merge'

workflow COMBINE {

    take:
    ch_h5ad  // queue channel: [ val(meta), path(h5ad) ]
    ch_base  // value channel: [ val(meta), path(h5ad) ]

    main:

    ch_versions      = channel.empty()
    ch_obs           = channel.empty()
    ch_var           = channel.empty()
    ch_obsm          = channel.empty()
    ch_celltypes     = channel.empty()

    ADATA_MERGE(
        ch_h5ad.map { _meta, h5ad -> [[id: "merged"], h5ad] }.groupTuple(),
        ch_base,
    )
    ch_var = ch_var.mix(ADATA_MERGE.out.intersect_genes)
    ch_outer = ADATA_MERGE.out.outer
    ch_inner = ADATA_MERGE.out.inner
    ch_versions = ch_versions.mix(ADATA_MERGE.out.versions)
    ch_celltypes = ADATA_MERGE.out.celltypes

    INTEGRATE(
        ADATA_MERGE.out.integrate,
        params.base_adata != null,
        params.integration_hvgs,
        params.integration_excluded_genes ? file(params.integration_excluded_genes) : [],
        params.integration_methods.split(',').collect { it -> it.trim().toLowerCase() },
        params.scvi_model,
        params.scanvi_model,
        params.scvi_categorical_covariates,
        params.scvi_continuous_covariates,
        params.scimilarity_model
    )
    ch_versions      = ch_versions.mix(INTEGRATE.out.versions)
    ch_var           = ch_var.mix(INTEGRATE.out.var)

    if (params.base_adata) {
        ADATA_MERGEEMBEDDINGS(
            INTEGRATE.out.integrations
            .combine(
                ch_base.map{ _meta, base -> base }
            ).combine(
                ADATA_MERGE.out.inner.map{ _meta, inner -> inner }
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
    h5ad             = ch_outer         // channel: [ merged, h5ad ]
    h5ad_inner       = ch_inner         // channel: [ merged, h5ad ]
    celltypes        = ch_celltypes     // channel: [ csv ]
    integrations     = ch_integrations  // channel: [ integration, h5ad ]
    var              = ch_var           // channel: [ pkl ]
    obs              = ch_obs           // channel: [ pkl ]
    obsm             = ch_obsm          // channel: [ pkl ]
    versions         = ch_versions      // channel: [ versions.yml ]
}
