include { SCANPY_HVGS        } from '../../../modules/local/scanpy/hvgs'
include { SCANPY_FILTER      } from '../../../modules/local/scanpy/filter'
include { SCVITOOLS_SCVI     } from '../../../modules/local/scvitools/scvi'
include { SCVITOOLS_SCANVI   } from '../../../modules/local/scvitools/scanvi'
include { SYMPHONY_HARMONYINTEGRATE } from '../../../modules/local/symphony/harmonyintegrate'
include { SYMPHONY_MAPEMBEDDING      } from '../../../modules/local/symphony/mapembedding'
include { SCANPY_BBKNN       } from '../../../modules/local/scanpy/bbknn'
include { SCANPY_COMBAT      } from '../../../modules/local/scanpy/combat'
include { SCANPY_PCA         } from '../../../modules/local/scanpy/pca'
include { SCARCHES_EXPIMAP   } from '../../../modules/local/scarches/expimap'
include { SEURAT_INTEGRATION } from '../../../modules/local/seurat/integration'
include { ADATA_READRDS      } from '../../../modules/local/adata/readrds'
include { SCIMILARITY        } from '../scimilarity'

def integrationMeta(meta, method) {
    def subset_suffix = meta.subset ? "-${meta.subset}" : ""
    meta + [
        id: "${method}${subset_suffix}",
        integration: method,
    ]
}

def integrationKey(meta) {
    meta.subset ?: 'all'
}

workflow INTEGRATE {
    take:
    ch_h5ad                     // channel: [ merged, h5ad ]
    is_extension                // boolean
    n_hvgs                      // integer
    excluded_genes              // path
    methods                     // list of string
    scvi_model                  // path
    scanvi_model                // path
    scvi_categorical_covariates // list of string
    scvi_continuous_covariates  // list of string
    scimilarity_model           // path
    symphony_reference           // path
    expimap_gmt                 // path
    condition_col               // string
    batch_col                   // string
    scanvi_label_col            // string
    scanvi_unlabeled_category   // string

    main:
    ch_versions = channel.empty()
    ch_obs = channel.empty()
    ch_var = channel.empty()
    ch_obsm = channel.empty()
    ch_integrations = channel.empty()

    // If a reference model is provided, only the genes in the reference model are used
    // Otherwise, we would intersect the HVGs, which is not what we want
    if (!is_extension) {
        SCANPY_HVGS (
            ch_h5ad,
            n_hvgs,
            excluded_genes
        )
        ch_versions = ch_versions.mix(SCANPY_HVGS.out.versions)
        ch_h5ad_hvg = SCANPY_HVGS.out.h5ad

        // See issue 215
        // ch_var = ch_var.mix(SCANPY_HVGS.out.var)

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
        ch_versions = ch_versions.mix(SCANPY_FILTER.out.versions)
    }
    else {
        ch_h5ad_hvg = ch_h5ad
    }

    if (methods.contains('seurat')) {
        SEURAT_INTEGRATION (
            ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'seurat'), h5ad] },
            batch_col
        )
        ch_versions = ch_versions.mix(SEURAT_INTEGRATION.out.versions)
        ch_integrations = ch_integrations.mix(SEURAT_INTEGRATION.out.h5ad)
    }

    if (methods.contains('scvi')) {
        SCVITOOLS_SCVI (
            (scvi_model ? ch_h5ad : ch_h5ad_hvg)
                .map { meta, h5ad -> [integrationMeta(meta, 'scvi'), h5ad] },
            scvi_model
                ? channel.value([[id: 'scvi'], scvi_model])
                : [[], []],
            batch_col,
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_versions = ch_versions.mix(SCVITOOLS_SCVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCVI.out.h5ad)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCVI.out.obsm)
    }

    if (methods.contains('scanvi')) {
        ch_scanvi_h5ad = (scvi_model ? ch_h5ad : ch_h5ad_hvg)
            .map { meta, h5ad -> [integrationMeta(meta, 'scanvi'), h5ad] }

        if (!scanvi_model && methods.contains('scvi')) {
            ch_scanvi_reference_model = ch_scanvi_h5ad
                .map { meta, h5ad -> [integrationKey(meta), meta, h5ad] }
                .join(
                    SCVITOOLS_SCVI.out.model
                        .map { meta, model -> [integrationKey(meta), meta, model] }
                )
                .multiMap { _key, meta, h5ad, meta2, model ->
                    h5ad: [meta, h5ad]
                    reference_model: [meta2, model]
                }
        }

        SCVITOOLS_SCANVI (
            !scanvi_model && methods.contains('scvi')
                ? ch_scanvi_reference_model.h5ad
                : ch_scanvi_h5ad,
            scanvi_model
                ? channel.value([[id: 'scanvi'], scanvi_model])
                : methods.contains('scvi')
                    ? ch_scanvi_reference_model.reference_model
                    : [[], []],
            [scanvi_label_col, scanvi_unlabeled_category],
            batch_col,
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_versions = ch_versions.mix(SCVITOOLS_SCANVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCANVI.out.h5ad)
        ch_obs = ch_obs.mix(SCVITOOLS_SCANVI.out.obs)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCANVI.out.obsm)
    }

    if (methods.contains('symphony')) {
        if (symphony_reference) {
            SYMPHONY_MAPEMBEDDING (
                ch_h5ad.map { meta, h5ad -> [integrationMeta(meta, 'symphony'), h5ad] },
                channel.value([[id: 'symphony'], symphony_reference]),
                batch_col,
                "X"
            )
            ch_versions = ch_versions.mix(SYMPHONY_MAPEMBEDDING.out.versions)
            ch_integrations = ch_integrations.mix(SYMPHONY_MAPEMBEDDING.out.h5ad)
            ch_obsm = ch_obsm.mix(SYMPHONY_MAPEMBEDDING.out.obsm)
        }
        else {
            SYMPHONY_HARMONYINTEGRATE (
                ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'symphony'), h5ad] },
                batch_col,
                "X"
            )
            ch_versions = ch_versions.mix(SYMPHONY_HARMONYINTEGRATE.out.versions)
            ch_integrations = ch_integrations.mix(SYMPHONY_HARMONYINTEGRATE.out.h5ad)
            ch_obsm = ch_obsm.mix(SYMPHONY_HARMONYINTEGRATE.out.obsm)
        }
    }

    if (methods.contains('bbknn')) {
        SCANPY_BBKNN (
            ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'bbknn'), h5ad] },
            batch_col
        )
        ch_versions = ch_versions.mix(SCANPY_BBKNN.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_BBKNN.out.h5ad)
    }

    if (methods.contains('combat')) {
        SCANPY_COMBAT (
            ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'combat'), h5ad] },
            batch_col
        )
        ch_versions = ch_versions.mix(SCANPY_COMBAT.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_COMBAT.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_COMBAT.out.obsm)
    }

    if (methods.contains('pca')) {
        SCANPY_PCA (
            ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'pca'), h5ad] },
            "X_emb"
        )
        ch_versions = ch_versions.mix(SCANPY_PCA.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_PCA.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_PCA.out.obsm)
    }

    if (methods.contains('expimap')) {
        SCARCHES_EXPIMAP (
            ch_h5ad_hvg.map { meta, h5ad -> [integrationMeta(meta, 'expimap'), h5ad] },
            expimap_gmt
            ? channel.value([[id: 'expimap'], file(expimap_gmt, checkIfExists: true)])
            : channel.value([[id: 'expimap'], file("${projectDir}/assets/databases/expimap/pathways.gmt", checkIfExists: true)]),
            condition_col,
            "X"
        )
        ch_versions = ch_versions.mix(SCARCHES_EXPIMAP.out.versions)
        ch_integrations = ch_integrations.mix(SCARCHES_EXPIMAP.out.h5ad)
        ch_obsm = ch_obsm.mix(SCARCHES_EXPIMAP.out.obsm)
    }

    if (methods.contains('scimilarity')) {
        SCIMILARITY (
            ch_h5ad.map { meta, h5ad -> [integrationMeta(meta, 'scimilarity'), h5ad] },
            scimilarity_model,
        )
        ch_versions = ch_versions.mix(SCIMILARITY.out.versions)
        ch_integrations = ch_integrations.mix(SCIMILARITY.out.integrations)
        ch_obs = ch_obs.mix(SCIMILARITY.out.obs)
        ch_obsm = ch_obsm.mix(SCIMILARITY.out.obsm)
    }

    emit:
    integrations = ch_integrations // channel: [ integration, h5ad ]
    obs          = ch_obs // channel: [ pkl ]
    var          = ch_var // channel: [ pkl ]
    obsm         = ch_obsm // channel: [ pkl ]
    versions     = ch_versions // channel: [ versions.yml ]
}
