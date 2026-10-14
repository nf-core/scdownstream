include { SCANPY_HVGS } from '../../../modules/local/scanpy/hvgs'
include { SCANPY_PEARSONRESIDUALS_HVGS } from '../../../modules/local/scanpy/pearsonresidualshvgs'
include { SCRY_DEVIANCE } from '../../../modules/local/scry/deviance'
include { SCANPY_FILTER } from '../../../modules/local/scanpy/filter'

workflow FEATURE_SELECTION {
    take:
    ch_h5ad // channel: [ meta, h5ad|zarr ]
    feature_selection //   value: string
    n_features //   value: integer
    excluded_genes //    path: file or []
    normalization_method //   value: string
    h5ad_required //   value: boolean

    main:
    ch_var = channel.empty()

    if (feature_selection == 'hvgs') {
        SCANPY_HVGS(
            ch_h5ad,
            n_features,
            excluded_genes,
            normalization_method,
            false,
            ch_h5ad.map { meta, _h5ad -> meta.batch_col ?: '' },
            'zarr',
        )
        ch_h5ad = SCANPY_HVGS.out.anndata
        ch_var = ch_var.mix(SCANPY_HVGS.out.var)
    }
    else if (feature_selection == 'deviance') {
        SCRY_DEVIANCE(
            ch_h5ad,
            n_features,
            excluded_genes,
            ch_h5ad.map { _meta, _h5ad -> _meta.batch_col ?: '' },
        )
        ch_h5ad = SCRY_DEVIANCE.out.h5ad
        ch_var = ch_var.mix(SCRY_DEVIANCE.out.var)
    }
    else if (feature_selection == 'pearson_residuals_hvgs') {
        SCANPY_PEARSONRESIDUALS_HVGS(
            ch_h5ad,
            n_features,
            excluded_genes,
            ch_h5ad.map { meta, _h5ad -> meta.batch_col ?: '' },
            ch_h5ad.map { meta, _h5ad -> meta.counts_layer ?: 'X' },
            'zarr',
        )
        ch_h5ad = SCANPY_PEARSONRESIDUALS_HVGS.out.anndata
        ch_var = ch_var.mix(SCANPY_PEARSONRESIDUALS_HVGS.out.var)
    }
    else if (feature_selection == 'none') {
        ch_h5ad = ch_h5ad
    }
    else {
        error("Unknown feature_selection: ${feature_selection}")
    }

    SCANPY_FILTER(
        ch_h5ad,
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
        [],
        false,
        h5ad_required ? 'h5ad' : 'zarr',
    )
    ch_h5ad = SCANPY_FILTER.out.anndata

    emit:
    anndata = ch_h5ad
    var = ch_var
}
