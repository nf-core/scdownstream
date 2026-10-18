include { CELDA_DECONTX               } from '../../../modules/local/celda/decontx'
include { CELLBENDER_REMOVEBACKGROUND } from '../../../modules/nf-core/cellbender/removebackground'
include { CELLBENDER_MERGE            } from '../../../modules/nf-core/cellbender/merge'
include { SOUPX                       } from '../../../modules/local/soupx'
include { SCVITOOLS_SCAR              } from '../../../modules/nf-core/scvitools/scar'

def resolve_ambient_expression(meta, method, global_flag) {
    def raw_layer = meta.counts_layer ?: 'X'
    def applied_method = meta.ambient_correction == false || method == 'none' ? 'none' : method
    def canonical_layer = 'ambient_corrected_' + applied_method
    def correction_layer = applied_method == 'none' ? null : canonical_layer
    if (correction_layer == raw_layer) {
        correction_layer = canonical_layer + '_internal'
    }
    def select_corrected = meta.ambient_corrected_integration == null
        ? global_flag == true
        : meta.ambient_corrected_integration == true
    select_corrected = select_corrected && applied_method != 'none'

    return meta + [
        raw_counts_layer: raw_layer,
        counts_layer: select_corrected ? correction_layer : raw_layer,
        ambient_correction_layer: correction_layer,
        ambient_correction_method: applied_method,
        analysis_counts_source: select_corrected ? canonical_layer : 'raw_counts'
    ]
}

workflow AMBIENT_CORRECTION {
    take:
    ch_pairing                    // channel: [ meta, h5ad, h5ad ]
    method                        // value: string
    ambient_corrected_integration // value: boolean
    scvi_max_epochs               // value: integer

    main:

    ch_do_ambient_correction = ch_pairing
        .branch {
            meta, _filtered, _unfiltered ->
            no: meta.ambient_correction == false
            yes: true
        }

    ch_multi = ch_do_ambient_correction.yes
        .multiMap {
            meta, filtered, unfiltered ->
            input: [meta, filtered, unfiltered]
            batch_col: meta.batch_col ?: []
            input_layer: meta.counts_layer ?: "X"
            output_layer: resolve_ambient_expression(meta, method, ambient_corrected_integration).ambient_correction_layer
        }

    ch_h5ad = ch_do_ambient_correction.no
        .map { meta, filtered, _unfiltered -> [meta, filtered] }

    if (method == 'none') {
        log.info("AMBIENT_CORRECTION: Not performed since 'none' selected.")
        ch_h5ad = ch_h5ad
            .mix(ch_multi.input
                .map { meta, filtered, _unfiltered -> [meta, filtered] })
    }
    else if (method == 'decontx') {
        CELDA_DECONTX (
            ch_multi.input,
            ch_multi.batch_col,
            ch_multi.input_layer,
            ch_multi.output_layer
        )
        ch_h5ad = ch_h5ad.mix(CELDA_DECONTX.out.h5ad)
    }
    else if (method == 'cellbender') {
        CELLBENDER_REMOVEBACKGROUND (
            ch_multi.input
                .map { meta, _filtered, unfiltered -> [meta, unfiltered] })

        CELLBENDER_MERGE (
            ch_multi.input
                .map { meta, filtered, raw -> [meta.id, meta, filtered, raw] }
                .join(CELLBENDER_REMOVEBACKGROUND.out.h5
                    .map {
                        meta, h5 ->
                        [meta.id, h5]
                    }, by: 0, failOnMismatch: true)
                .map {
                    _id, meta, filtered, raw, h5 ->
                    [meta, filtered, raw, h5]
                },
            ch_multi.output_layer
        )
        ch_h5ad = ch_h5ad.mix(CELLBENDER_MERGE.out.h5ad)
    }
    else if (method == 'soupx') {
        SOUPX (
            ch_multi.input,
            0.8,
            ch_multi.input_layer,
            ch_multi.output_layer
        )
        ch_h5ad = ch_h5ad.mix(SOUPX.out.h5ad)
    }
    else if (method == 'scar') {
        SCVITOOLS_SCAR (
            ch_multi.input,
            ch_multi.input_layer,
            ch_multi.output_layer,
            scvi_max_epochs ?: [],
            []
        )
        ch_h5ad = ch_h5ad.mix(SCVITOOLS_SCAR.out.h5ad)
    }
    else {
        error("AMBIENT_CORRECTION: Unexpected method for ambient RNA correction: '${method}'.")
    }

    ch_h5ad = ch_h5ad.map { meta, h5ad ->
        [resolve_ambient_expression(meta, method, ambient_corrected_integration), h5ad]
    }

    emit:
    h5ad     = ch_h5ad     // channel: [ meta, h5ad ]
}
