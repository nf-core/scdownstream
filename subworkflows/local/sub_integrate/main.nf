include { ADATA_SPLITCOL as SPLITCOL } from '../../../modules/local/adata/splitcol'
include { INTEGRATE                  } from '../integrate'

workflow SUB_INTEGRATE {
    take:
    ch_h5ad                     // channel: [ val(meta), path(h5ad) ]
    split_col                   //   value: string
    n_hvgs                      //   value: integer
    excluded_genes              //    path: file or []
    methods                     //   value: list of string
    scvi_model                  //   value: string
    scanvi_model                //   value: string
    scvi_categorical_covariates //   value: string
    scvi_continuous_covariates  //   value: string
    scimilarity_model           //   value: string
    symphony_reference          //    path: file or null
    expimap_gmt                 //   value: string
    condition_col               //   value: string
    label_whitelist             //   value: string or null

    main:
    def normalized_whitelist = label_whitelist
        ? label_whitelist.split(',')*.trim().findAll { label -> label }.collect { label -> label.replace(' ', '_') }
        : []

    SPLITCOL (
        ch_h5ad,
        split_col
    )

    ch_h5ad_split = SPLITCOL.out.h5ad
        .transpose()
        .map { meta, h5ad ->
            def subset = h5ad.simpleName
            [
                meta + [
                    id: subset,
                    subset: subset,
                ],
                h5ad,
            ]
        }

    if (normalized_whitelist) {
        ch_h5ad_split = ch_h5ad_split
            .filter { meta, _h5ad -> meta.subset in normalized_whitelist }
            .ifEmpty {
                error("integrate_per_label_whitelist: none of the requested labels matched any group in '${split_col}': ${normalized_whitelist.join(', ')}")
            }
    }

    INTEGRATE (
        ch_h5ad_split,
        false,
        n_hvgs,
        excluded_genes,
        methods,
        scvi_model,
        scanvi_model,
        scvi_categorical_covariates,
        scvi_continuous_covariates,
        scimilarity_model,
        symphony_reference,
        expimap_gmt,
        condition_col
    )

    ch_integrations = INTEGRATE.out.integrations
        .map { meta, h5ad ->
            [meta + [id: "${meta.integration}-${meta.subset}"], h5ad]
        }

    emit:
    integrations = ch_integrations            // channel: [ meta, h5ad ]
    obs          = INTEGRATE.out.obs          // channel: [ pkl ]
    var          = INTEGRATE.out.var          // channel: [ pkl ]
    obsm         = INTEGRATE.out.obsm         // channel: [ pkl ]
}
