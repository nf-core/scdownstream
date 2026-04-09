include { SCANPY_HVGS        } from '../../../modules/local/scanpy/hvgs'
include { SCANPY_FILTER      } from '../../../modules/local/scanpy/filter'
include { SCVITOOLS_SCVI     } from '../../../modules/local/scvitools/scvi'
include { SCVITOOLS_SCANVI   } from '../../../modules/local/scvitools/scanvi'
include { SCANPY_HARMONY     } from '../../../modules/local/scanpy/harmony'
include { SCANPY_BBKNN       } from '../../../modules/local/scanpy/bbknn'
include { SCANPY_COMBAT      } from '../../../modules/local/scanpy/combat'
include { SEURAT_INTEGRATION } from '../../../modules/local/seurat/integration'
include { ADATA_READRDS      } from '../../../modules/local/adata/readrds'
include { SCIMILARITY        } from '../scimilarity'

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
            ch_h5ad_hvg.map { _meta, h5ad -> [[id: 'seurat'], h5ad] }, "batch"
        )
        ch_versions = ch_versions.mix(SEURAT_INTEGRATION.out.versions)
        ch_integrations = ch_integrations.mix(SEURAT_INTEGRATION.out.h5ad)
    }

    if (methods.contains('scvi')) {
        SCVITOOLS_SCVI (
            (scvi_model ? ch_h5ad : ch_h5ad_hvg)
                .map { _meta, h5ad -> [[id: 'scvi'], h5ad] },
            scvi_model
                ? channel.value([[id: 'scvi'], scvi_model])
                : [[], []],
            "batch",
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_versions = ch_versions.mix(SCVITOOLS_SCVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCVI.out.h5ad)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCVI.out.obsm)
    }

    if (methods.contains('scanvi')) {
        SCVITOOLS_SCANVI (
            (scvi_model ? ch_h5ad : ch_h5ad_hvg)
                .map { _meta, h5ad -> [[id: 'scanvi'], h5ad] },
            scanvi_model
                ? channel.value([[id: 'scanvi'], scanvi_model])
                : methods.contains('scvi')
                    ? SCVITOOLS_SCVI.out.model
                    : [[], []],
            ["label", "Unknown"],
            "batch",
            scvi_categorical_covariates,
            scvi_continuous_covariates,
        )
        ch_versions = ch_versions.mix(SCVITOOLS_SCANVI.out.versions)
        ch_integrations = ch_integrations.mix(SCVITOOLS_SCANVI.out.h5ad)
        ch_obs = ch_obs.mix(SCVITOOLS_SCANVI.out.obs)
        ch_obsm = ch_obsm.mix(SCVITOOLS_SCANVI.out.obsm)
    }

    if (methods.contains('harmony')) {
        SCANPY_HARMONY (
            ch_h5ad_hvg.map { _meta, h5ad -> [[id: 'harmony'], h5ad] },
            "batch",
            "X"
        )
        ch_versions = ch_versions.mix(SCANPY_HARMONY.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_HARMONY.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_HARMONY.out.obsm)
    }

    if (methods.contains('bbknn')) {
        SCANPY_BBKNN (
            ch_h5ad_hvg.map { _meta, h5ad -> [[id: 'bbknn'], h5ad] },
            "batch"
        )
        ch_versions = ch_versions.mix(SCANPY_BBKNN.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_BBKNN.out.h5ad)
    }

    if (methods.contains('combat')) {
        SCANPY_COMBAT (
            ch_h5ad_hvg.map { _meta, h5ad -> [[id: 'combat'], h5ad] },
            "batch"
        )
        ch_versions = ch_versions.mix(SCANPY_COMBAT.out.versions)
        ch_integrations = ch_integrations.mix(SCANPY_COMBAT.out.h5ad)
        ch_obsm = ch_obsm.mix(SCANPY_COMBAT.out.obsm)
    }

    if (methods.contains('scimilarity')) {
        SCIMILARITY (
            ch_h5ad.map { _meta, h5ad -> [[id: 'scimilarity'], h5ad] },
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
