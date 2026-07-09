include { SCANPY_LOG1PNORM        } from '../../../modules/local/scanpy/log1pnorm'
include { SCANPY_PEARSONRESIDUALS } from '../../../modules/local/scanpy/pearsonresiduals'
include { SCRAN_NORMALIZATION     } from '../../../modules/local/scran/normalization'

workflow NORMALIZATION {
    take:
    ch_h5ad                 // channel: [ meta, h5ad ]
    normalization_methods   //   value: list of string
    transformed_layer       //   value: string

    main:
    ch_layers = channel.empty()
    ch_working = ch_h5ad

    if (normalization_methods.contains('log1p')) {
        SCANPY_LOG1PNORM(ch_h5ad)
        ch_layers = ch_layers.mix(SCANPY_LOG1PNORM.out.layers)
        if (transformed_layer == 'log1p_norm') {
            ch_working = SCANPY_LOG1PNORM.out.h5ad
        }
    }

    if (normalization_methods.contains('scran')) {
        SCRAN_NORMALIZATION(ch_h5ad)
        ch_layers = ch_layers.mix(SCRAN_NORMALIZATION.out.layers)
        if (transformed_layer == 'scran') {
            ch_working = SCRAN_NORMALIZATION.out.h5ad
        }
    }

    if (normalization_methods.contains('pearson_residuals')) {
        SCANPY_PEARSONRESIDUALS(ch_h5ad)
        ch_layers = ch_layers.mix(SCANPY_PEARSONRESIDUALS.out.layers)
        if (transformed_layer == 'pearson_residuals') {
            ch_working = SCANPY_PEARSONRESIDUALS.out.h5ad
        }
    }

    emit:
    h5ad   = ch_working
    layers = ch_layers
}
