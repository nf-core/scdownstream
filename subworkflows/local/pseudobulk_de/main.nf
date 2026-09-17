include { PYDESEQ2_DIFFERENTIAL } from '../../../modules/local/pydeseq2/differential'
include { EDGEPYTHON_DIFFERENTIAL } from '../../../modules/local/edgepython/differential'

workflow PSEUDOBULK_DE {
    take:
    ch_h5ad_pseudobulk     // channel: [ meta, h5ad ] with meta.de_methods_resolved
    reference_condition    //   value: string
    interesting_genes      //   value: string (path) or []

    main:
    ch_results       = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_pydeseq2 = ch_h5ad_pseudobulk
        .filter { meta, _h5ad -> 'pydeseq2' in meta.de_methods_resolved }

    ch_edgepython = ch_h5ad_pseudobulk
        .filter { meta, _h5ad -> 'edgepython' in meta.de_methods_resolved }

    PYDESEQ2_DIFFERENTIAL(
        ch_pydeseq2,
        '',
        reference_condition ?: '',
        interesting_genes ?: [],
    )
    ch_results       = ch_results.mix(PYDESEQ2_DIFFERENTIAL.out.results)
    ch_multiqc_files = ch_multiqc_files.mix(PYDESEQ2_DIFFERENTIAL.out.multiqc_files)

    EDGEPYTHON_DIFFERENTIAL(
        ch_edgepython,
        reference_condition ?: '',
        interesting_genes ?: [],
    )
    ch_results       = ch_results.mix(EDGEPYTHON_DIFFERENTIAL.out.results)
    ch_multiqc_files = ch_multiqc_files.mix(EDGEPYTHON_DIFFERENTIAL.out.multiqc_files)

    emit:
    results       = ch_results
    multiqc_files = ch_multiqc_files
}
