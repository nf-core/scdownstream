include { ADATA_SPLITCOL as SPLITCOL } from '../../../modules/local/adata/splitcol'
include { INTEGRATE                  } from '../integrate'
include { anndata                      } from 'plugin/nf-anndata'

workflow SUB_INTEGRATE {
    take:
    ch_h5ad                     // channel: [ val(meta), path(h5ad) ]
    split_col                   //   value: string
    feature_selection           //   value: string
    n_features                  //   value: integer
    excluded_genes              //    path: file or []
    normalization_method        //   value: string
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
    def normalized_methods = methods*.trim()*.toLowerCase()

    def provided_reference_models = []
    if (scvi_model && normalized_methods.contains('scvi')) {
        provided_reference_models << 'scvi_model'
    }
    if (scanvi_model && normalized_methods.contains('scanvi')) {
        provided_reference_models << 'scanvi_model'
    }
    if (scimilarity_model && normalized_methods.contains('scimilarity')) {
        provided_reference_models << 'scimilarity_model'
    }
    if (symphony_reference && normalized_methods.contains('symphony')) {
        provided_reference_models << 'symphony_reference'
    }

    if (provided_reference_models) {
        log.warn """\
            Per-label integration (integrate_per_label) was enabled, but reference model parameter(s) [${provided_reference_models.join(', ')}] were also provided.
            Query cells will be mapped into the pre-trained reference latent space instead of training a separate integration model per label group.
            The resulting embeddings may therefore be identical across label groups. \
            """.stripIndent()
    }

    def normalized_whitelist = label_whitelist
        ? label_whitelist.split(',')*.trim().findAll { label -> label }.collect { label -> label.replace(' ', '_') }
        : []

    if (normalized_whitelist) {
        ch_h5ad = ch_h5ad.map { meta, h5ad ->
            def ad = anndata(h5ad)
            if (!(split_col in ad.obs.colnames)) {
                error("integrate_per_label_whitelist: column '${split_col}' not found in adata")
            }
            def available_groups = ad.obs[split_col].unique().collect { value -> value.toString().replace(' ', '_') }
            if (!normalized_whitelist.any { label -> label in available_groups }) {
                error("integrate_per_label_whitelist: none of the requested labels matched any group in '${split_col}': ${normalized_whitelist.join(', ')}")
            }
            [meta, h5ad]
        }
    }

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
    }

    INTEGRATE (
        ch_h5ad_split,
        false,
        feature_selection,
        n_features,
        excluded_genes,
        normalization_method,
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
    layers       = INTEGRATE.out.layers       // channel: [ *.npy ]
}
