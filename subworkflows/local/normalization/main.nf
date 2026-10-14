include { SCANPY_LOG1PNORM } from '../../../modules/local/scanpy/log1pnorm'
include { SCRAN_NORMALIZATION } from '../../../modules/local/scran/normalization'

workflow NORMALIZATION {
    take:
    ch_h5ad // channel: [ meta, h5ad|zarr ]
    normalization_method //   value: string

    main:
    if (normalization_method == 'log1p') {
        SCANPY_LOG1PNORM(ch_h5ad, 'zarr')
        ch_h5ad = SCANPY_LOG1PNORM.out.anndata
    }
    else if (normalization_method == 'scran') {
        SCRAN_NORMALIZATION(ch_h5ad)
        ch_h5ad = SCRAN_NORMALIZATION.out.h5ad
    }
    else {
        error("Unknown normalization_method: ${normalization_method}. Supported values: log1p, scran")
    }

    emit:
    anndata = ch_h5ad // channel: [ meta, h5ad|zarr ]
}
