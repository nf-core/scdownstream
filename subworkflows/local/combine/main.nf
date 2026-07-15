include { INTEGRATE                              } from '../integrate'
include { NORMALIZATION                          } from '../normalization'
include { ADATA_MERGEEMBEDDINGS                  } from '../../../modules/local/adata/mergeembeddings'
include { ADATA_MERGE                            } from '../../../modules/local/adata/merge'
include { ADATA_SUBSETVAR as SUBSETVAR_INNER     } from '../../../modules/local/adata/subsetvar'
include { ADATA_SUBSETVAR as SUBSETVAR_INTEGRATE } from '../../../modules/local/adata/subsetvar'
include { SCIBMETRICS_BENCHMARK                  } from '../../../modules/local/scibmetrics/benchmark'

workflow COMBINE {

    take:
    ch_h5ad                     // channel: [ val(meta), path(h5ad) ]
    ch_base                     // channel: [ val(meta), path(h5ad) ]
    is_extension                //   value: boolean
    feature_selection           //   value: string
    integration_n_features      //   value: integer
    integration_methods         //   value: string
    integration_excluded_genes  //   value: string
    integration_exclude_mt      //   value: boolean
    mito_genes                  //   value: string
    normalization_method        //   value: string
    scvi_model                  //   value: string
    scanvi_model                //   value: string
    scvi_categorical_covariates //   value: string
    scvi_continuous_covariates  //   value: string
    scimilarity_model           //   value: string
    symphony_reference          //   value: string
    expimap_gmt                 //   value: string
    condition_col               //   value: string
    scib                        //   value: boolean
    scib_max_cells              //   value: integer or null
    scib_subsample_strategy     //   value: string
    scib_subsample_seed         //   value: integer
    scib_metric_profile         //   value: string

    main:

    ch_multiqc_files = channel.empty()
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
    ch_outer = ADATA_MERGE.out.outer

    if (!is_extension) {
        NORMALIZATION(
            ch_outer,
            normalization_method,
        )
        ch_outer = NORMALIZATION.out.h5ad

        SUBSETVAR_INNER(
            ch_outer,
            'intersection',
        )
        ch_inner = SUBSETVAR_INNER.out.h5ad
        ch_integrate_input = ch_inner
    }
    else {
        SUBSETVAR_INNER(
            ch_outer,
            'intersection',
        )
        ch_inner = SUBSETVAR_INNER.out.h5ad

        SUBSETVAR_INTEGRATE(
            ADATA_MERGE.out.integrate,
            'intersection',
        )
        ch_integrate_input = SUBSETVAR_INTEGRATE.out.h5ad
    }

    INTEGRATE(
        ch_integrate_input,
        is_extension,
        feature_selection,
        integration_n_features,
        integration_excluded_genes ? file(integration_excluded_genes) : [],
        integration_exclude_mt,
        mito_genes ? file(mito_genes) : [],
        normalization_method,
        integration_methods
            .split(',')
            .collect { it -> it.trim().toLowerCase() },
        scvi_model,
        scanvi_model,
        scvi_categorical_covariates,
        scvi_continuous_covariates,
        scimilarity_model,
        symphony_reference,
        expimap_gmt,
        condition_col
    )
    ch_var = INTEGRATE.out.var

    if (is_extension) {
        ADATA_MERGEEMBEDDINGS(
            INTEGRATE.out.integrations
                .map { meta, integrated -> [meta, meta.integration ?: meta.id, integrated] }
                .combine(
                    ch_base.map { _meta, base -> base }
                ).combine(
                    ch_inner.map { _meta, inner -> inner }
                )
        )
        ch_integrations  = ADATA_MERGEEMBEDDINGS.out.h5ad
        ch_obs           = ch_obs.mix(ADATA_MERGEEMBEDDINGS.out.obs)
        ch_obsm          = ch_obsm.mix(ADATA_MERGEEMBEDDINGS.out.obsm)
    } else {
        ch_integrations  = INTEGRATE.out.integrations
        ch_obs           = ch_obs.mix(INTEGRATE.out.obs)
        ch_obsm          = ch_obsm.mix(INTEGRATE.out.obsm)
    }

    ch_integrations = ch_integrations
        .map { meta, file -> [meta + [id: meta.integration], file] }

    if (scib) {
        SCIBMETRICS_BENCHMARK (
            ch_integrations
                // BBKNN corrects the neighborhood graph and does not produce a dense embedding
                // Thus, it is not compatible with scib-metrics
                .filter { meta, _h5ad -> meta.integration != 'bbknn' },
            scib_max_cells ?: 0,
            scib_subsample_strategy,
            scib_subsample_seed,
            scib_metric_profile,
        )
        ch_multiqc_files = ch_multiqc_files.mix(SCIBMETRICS_BENCHMARK.out.multiqc_files)
    }

    emit:
    h5ad             = ch_outer         // channel: [ merged, h5ad ]
    h5ad_inner       = ch_inner         // channel: [ merged, h5ad ]
    integrations     = ch_integrations  // channel: [ integration, h5ad ]
    var              = ch_var           // channel: [ pkl ]
    obs              = ch_obs           // channel: [ pkl ]
    obsm             = ch_obsm          // channel: [ pkl ]
    multiqc_files    = ch_multiqc_files // channel: [ *_mqc.json ]
}
