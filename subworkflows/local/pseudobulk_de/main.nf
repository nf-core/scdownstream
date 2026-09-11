include { PYDESEQ2_DIFFERENTIAL   } from '../../../modules/local/pydeseq2/differential'
include { EDGEPYTHON_DIFFERENTIAL } from '../../../modules/local/edgepython/differential'

workflow PSEUDOBULK_DE {
    take:
    ch_h5ad_pseudobulk  // channel: [ meta, h5ad ] with meta.de_methods_resolved
    reference_condition //   value: string

    main:
    ch_results = channel.empty()

    ch_pydeseq2 = ch_h5ad_pseudobulk
        .filter { meta, _h5ad -> 'pydeseq2' in meta.de_methods_resolved }
        .map { meta, h5ad -> [meta + [de_method: 'pydeseq2'], h5ad] }

    ch_edgepython = ch_h5ad_pseudobulk
        .filter { meta, _h5ad -> 'edgepython' in meta.de_methods_resolved }
        .map { meta, h5ad -> [meta + [de_method: 'edgepython'], h5ad] }

    PYDESEQ2_DIFFERENTIAL(
        ch_pydeseq2,
        '',
        reference_condition ?: '',
    )
    ch_results = ch_results.mix(PYDESEQ2_DIFFERENTIAL.out.results)

    EDGEPYTHON_DIFFERENTIAL(
        ch_edgepython,
        reference_condition ?: '',
    )
    ch_results = ch_results.mix(EDGEPYTHON_DIFFERENTIAL.out.results)

    emit:
    results = ch_results
}
