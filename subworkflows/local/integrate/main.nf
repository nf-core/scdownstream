include { SCANPY_HVGS        } from '../../../modules/local/scanpy/hvgs'
include { SCRY_DEVIANCE       } from '../../../modules/local/scry/deviance'
include { SCANPY_FILTER      } from '../../../modules/local/scanpy/filter'
include { NORMALIZATION      } from '../normalization'
include { SCVITOOLS_SCVI     } from '../../../modules/local/scvitools/scvi'
include { SCVITOOLS_SCANVI   } from '../../../modules/local/scvitools/scanvi'
include { SYMPHONY_HARMONYINTEGRATE } from '../../../modules/local/symphony/harmonyintegrate'
include { SYMPHONY_MAPEMBEDDING     } from '../../../modules/local/symphony/mapembedding'
include { SCANPY_BBKNN       } from '../../../modules/local/scanpy/bbknn'
include { SCANPY_COMBAT      } from '../../../modules/local/scanpy/combat'
include { SCANPY_PCA         } from '../../../modules/local/scanpy/pca'
include { SCARCHES_EXPIMAP   } from '../../../modules/local/scarches/expimap'
include { SEURAT_INTEGRATION } from '../../../modules/local/seurat/integration'
include { ADATA_READRDS      } from '../../../modules/local/adata/readrds'
include { SCIMILARITY        } from '../scimilarity'

workflow INTEGRATE {
    take:
    ch_h5ad                     // channel: [ merged, h5ad ]
    is_extension                // boolean
    feature_selection           // string: hvgs | deviance | none
    n_features                      // integer
    excluded_genes              // path
    normalization_methods       // list of string
    transformed_layer           // string
    methods                     // list of string
    scvi_model                  // path
    scanvi_model                // path
    scvi_categorical_covariates // list of string
    scvi_continuous_covariates  // list of string
    scimilarity_model           // path
    symphony_reference          // path
    expimap_gmt                 // path
    condition_col               // string

    main:
    ch_obs = channel.empty()
    ch_var = channel.empty()
    ch_obsm = channel.empty()
    ch_layers = channel.empty()
    ch_integrations = channel.empty()

    // If a reference model is provided, only the genes in the reference model are used
    // Otherwise, we would intersect the HVGs, which is not what we want
    if (!is_extension) {
        if (normalization_methods) {
            NORMALIZATION(
                ch_h5ad,
                normalization_methods,
                transformed_layer ?: '',
            )
            ch_h5ad = NORMALIZATION.out.h5ad
            ch_layers = ch_layers.mix(NORMALIZATION.out.layers)
        }

        if (feature_selection == 'hvgs') {
            SCANPY_HVGS (
                ch_h5ad,
                n_features,
                excluded_genes,
                transformed_layer ?: '',
            )
            ch_h5ad_hvg = SCANPY_HVGS.out.h5ad
            ch_var = ch_var.mix(SCANPY_HVGS.out.var)
        }
        else if (feature_selection == 'deviance') {
            SCRY_DEVIANCE (
                ch_h5ad,
                n_features,
                excluded_genes,
                ch_h5ad.map { _meta, _h5ad -> _meta.batch_col ?: '' }
            )
            ch_h5ad_hvg = SCRY_DEVIANCE.out.h5ad
            ch_var = ch_var.mix(SCRY_DEVIANCE.out.var)
        }
        else if (feature_selection == 'none') {
            ch_h5ad_hvg = ch_h5ad
        }
        else {
            error("Unknown feature_selection: ${feature_selection}")
        }

        // Filter out empty cells from the AnnData object
        SCANPY_FILTER (
            ch_h5ad_hvg,
            "index",
            1,
            0,
            0,
            0,
            100,
            0,
            100,
            0,
            0,
            0,
            0,
            []
        )
        ch_h5ad_hvg = SCANPY_FILTER.out.h5ad
    }
    else {
        ch_h5ad_hvg = ch_h5ad
    }

    if (methods.contains('seurat')) {
        SEURAT_INTEGRATION (
            ch_h5ad_hvg.map { meta, h5ad ->
                [meta + [integration: 'seurat'], h5ad]
            },
            "batch"
        )
        ch_integrations = ch_integrations.mix(SEURAT_INTEGRATION.out.h5ad)
    }

    if (methods.contains('scvi')) {
        SCVITOOLS_SCVI (
            (scvi_model ? ch_h5ad : ch_h5ad_hvg)
                .map { meta, h5ad ->
                    [meta + [integration: 'scvi'], h5ad]
                },
            scvi_model
                ? channel.value([[id: 'scvi'], scvi_model])
                : [[], []],
            "batch",
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCVI.out.h5ad)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCVI.out.obsm)
    }

    if (methods.contains('scanvi')) {
        ch_scanvi_h5ad = (scvi_model ? ch_h5ad : ch_h5ad_hvg)
            .map { meta, h5ad ->
                [meta + [integration: 'scanvi'], h5ad]
            }

        if (!scanvi_model && methods.contains('scvi')) {
            ch_scanvi_reference_model = ch_scanvi_h5ad
                .map { meta, h5ad -> [meta.id, meta, h5ad] }
                .join(
                    SCVITOOLS_SCVI.out.model
                        .map { meta, model -> [meta.id, meta, model] }
                )
                .multiMap { _key, meta, h5ad, meta2, model ->
                    h5ad: [meta, h5ad]
                    reference_model: [meta2, model]
                    reference_model_type: meta2.integration ?: meta2.id ?: ''
                }
        }

        ch_scanvi_reference_model_type = !scanvi_model && methods.contains('scvi')
            ? ch_scanvi_reference_model.reference_model_type
            : scanvi_model
                ? ch_scanvi_h5ad.map { _meta, _h5ad -> 'scanvi' }
                : ch_scanvi_h5ad.map { _meta, _h5ad -> '' }

        SCVITOOLS_SCANVI (
            !scanvi_model && methods.contains('scvi')
                ? ch_scanvi_reference_model.h5ad
                : ch_scanvi_h5ad,
            scanvi_model
                ? channel.value([[id: 'scanvi'], scanvi_model])
                : methods.contains('scvi')
                    ? ch_scanvi_reference_model.reference_model
                    : [[], []],
            ch_scanvi_reference_model_type,
            ["label", "Unknown"],
            "batch",
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCANVI.out.h5ad)
        ch_obs = ch_obs.mix(SCVITOOLS_SCANVI.out.obs)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCANVI.out.obsm)
    }

    if (methods.contains('symphony')) {
        if (symphony_reference) {
            SYMPHONY_MAPEMBEDDING (
                ch_h5ad.map { meta, h5ad ->
                    [meta + [integration: 'symphony'], h5ad]
                },
                channel.value([[id: 'symphony'], symphony_reference]),
                "batch",
                "X"
            )
            ch_integrations = ch_integrations.mix(SYMPHONY_MAPEMBEDDING.out.h5ad)
            ch_obsm = ch_obsm.mix(SYMPHONY_MAPEMBEDDING.out.obsm)
        }
        else {
            SYMPHONY_HARMONYINTEGRATE (
                ch_h5ad_hvg.map { meta, h5ad ->
                    [meta + [integration: 'symphony'], h5ad]
                },
                "batch",
                "X"
            )
            ch_integrations = ch_integrations.mix(SYMPHONY_HARMONYINTEGRATE.out.h5ad)
            ch_obsm = ch_obsm.mix(SYMPHONY_HARMONYINTEGRATE.out.obsm)
        }
    }

    if (methods.contains('bbknn')) {
        SCANPY_BBKNN (
            ch_h5ad_hvg.map { meta, h5ad ->
                [meta + [integration: 'bbknn'], h5ad]
            },
            "batch"
        )
        ch_integrations = ch_integrations.mix(SCANPY_BBKNN.out.h5ad)
    }

    if (methods.contains('combat')) {
        SCANPY_COMBAT (
            ch_h5ad_hvg.map { meta, h5ad ->
                [meta + [integration: 'combat'], h5ad]
            },
            "batch"
        )
        ch_integrations = ch_integrations.mix(SCANPY_COMBAT.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_COMBAT.out.obsm)
    }

    if (methods.contains('pca')) {
        SCANPY_PCA (
            ch_h5ad_hvg.map { meta, h5ad ->
                [meta + [integration: 'pca'], h5ad]
            },
            "X_emb",
            transformed_layer ?: '',
        )
        ch_integrations = ch_integrations.mix(SCANPY_PCA.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_PCA.out.obsm)
    }

    if (methods.contains('expimap')) {
        SCARCHES_EXPIMAP (
            ch_h5ad_hvg.map { meta, h5ad ->
                [meta + [integration: 'expimap'], h5ad]
            },
            expimap_gmt
            ? channel.value([[id: 'expimap'], file(expimap_gmt, checkIfExists: true)])
            : channel.value([[id: 'expimap'], file("${projectDir}/assets/databases/expimap/pathways.gmt", checkIfExists: true)]),
            condition_col,
            "X"
        )
        ch_integrations = ch_integrations.mix(SCARCHES_EXPIMAP.out.h5ad)
        ch_obsm = ch_obsm.mix(SCARCHES_EXPIMAP.out.obsm)
    }

    if (methods.contains('scimilarity')) {
        SCIMILARITY (
            ch_h5ad.map { meta, h5ad ->
                [meta + [integration: 'scimilarity'], h5ad]
            },
            scimilarity_model,
        )
        ch_integrations = ch_integrations.mix(SCIMILARITY.out.integrations)
        ch_obs = ch_obs.mix(SCIMILARITY.out.obs)
        ch_obsm = ch_obsm.mix(SCIMILARITY.out.obsm)
    }

    emit:
    integrations = ch_integrations // channel: [ integration, h5ad ]
    obs          = ch_obs // channel: [ pkl ]
    var          = ch_var // channel: [ pkl ]
    obsm         = ch_obsm // channel: [ pkl ]
    layers       = ch_layers // channel: [ *.npy ]
}
