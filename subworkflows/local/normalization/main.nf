include { SCANPY_LOG1PNORM    } from '../../../modules/local/scanpy/log1pnorm'
include { SCRAN_NORMALIZATION } from '../../../modules/local/scran/normalization'

workflow NORMALIZATION {
    take:
    ch_h5ad               // channel: [ meta, h5ad ]
    normalization_method  //   value: string

    main:
    if (normalization_method == 'log1p') {
        SCANPY_LOG1PNORM(ch_h5ad)
        ch_h5ad = SCANPY_LOG1PNORM.out.h5ad
        ch_layers = SCANPY_LOG1PNORM.out.layers
    }
    else if (normalization_method == 'scran') {
        SCRAN_NORMALIZATION(ch_h5ad)
        ch_h5ad = SCRAN_NORMALIZATION.out.h5ad
        ch_layers = SCRAN_NORMALIZATION.out.layers
    }
    else {
        error("Unknown normalization_method: ${normalization_method}. Supported values: log1p, scran")
    }

    emit:
    h5ad   = ch_h5ad
    layers = ch_layers
}
