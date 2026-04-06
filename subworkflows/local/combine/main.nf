include { INTEGRATE             } from '../integrate'
include { ADATA_MERGEEMBEDDINGS } from '../../../modules/local/adata/mergeembeddings'
include { ADATA_MERGE           } from '../../../modules/local/adata/merge'

workflow COMBINE {

    take:
    ch_h5ad                     // channel: [ val(meta), path(h5ad) ]
    ch_base                     // channel: [ val(meta), path(h5ad) ]
    is_extension                //   value: boolean
    integration_hvgs            //   value: integer
    integration_methods         //   value: string
    integration_excluded_genes  //   value: string
    scvi_model                  //   value: string
    scanvi_model                //   value: string
    scvi_categorical_covariates //   value: string
    scvi_continuous_covariates  //   value: string
    scimilarity_model           //   value: string

    main:

    ch_versions      = channel.empty()
    ch_obs           = channel.empty()
    ch_var           = channel.empty()
    ch_obsm          = channel.empty()

    ADATA_MERGE(
        ch_h5ad
            .map { _meta, h5ad -> [[id: "merged"], h5ad] }
            .groupTuple()
            .map { meta, h5ads -> [meta, h5ads.sort { a, b -> a.name <=> b.name }] },
        ch_base,
    )
    ch_var = ch_var.mix(ADATA_MERGE.out.intersect_genes)
    ch_outer = ADATA_MERGE.out.outer
    ch_inner = ADATA_MERGE.out.inner
    ch_versions = ch_versions.mix(ADATA_MERGE.out.versions)

    INTEGRATE(
        ADATA_MERGE.out.integrate,
        is_extension,
        integration_hvgs,
        integration_excluded_genes ? file(integration_excluded_genes) : [],
        integration_methods
            .split(',')
            .collect { it -> it.trim().toLowerCase() },
        scvi_model,
        scanvi_model,
        scvi_categorical_covariates,
        scvi_continuous_covariates,
        scimilarity_model
    )
    ch_versions      = ch_versions.mix(INTEGRATE.out.versions)
    ch_var           = ch_var.mix(INTEGRATE.out.var)

    if (is_extension) {
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
    integrations     = ch_integrations  // channel: [ integration, h5ad ]
    var              = ch_var           // channel: [ pkl ]
    obs              = ch_obs           // channel: [ pkl ]
    obsm             = ch_obsm          // channel: [ pkl ]
    versions         = ch_versions      // channel: [ versions.yml ]
}
