include { ADATA_TORDS } from '../../modules/local/adata/tords'

workflow GRN {
    take:
    ch_h5ad // channel: [ val(meta), path(h5ad) ]

    main:
    ch_versions = Channel.empty()

    if (prams.grn_methods == 'none') {
        log.info "GRN: Not performed since 'none' selected."
    } else {
        ADATA_TORDS(ch_h5ad)
        ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

        METACELLS(ADATA_TORDS.out.rds)
        ch_versions = ch_versions.mix(METACELLS.out.versions)
        ch_metacells = METACELLS.out.metacells

        methods = params.grn_methods.split(',').collect{it.trim().toLowerCase()}
        ch_grn_results = Channel.empty()

        if (methods.contains('boostdiff')) {
            BOOSTDIFF(ch_metacells)
            ch_versions = ch_versions.mix(BOOSTDIFF.out.versions)
            ch_grn_results = ch_grn_results.mix(BOOSTDIFF.out.grn_results)
        }

        if (methods.contains('zscores')) {
            ZSCORES(ch_metacells)
            ch_versions = ch_versions.mix(ZSCORES.out.versions)
            ch_grn_results = ch_grn_results.mix(ZSCORES.out.grn_results)
        }

        if (methods.contains('diffcoex')) {
            DIFFCOEX(ch_metacells)
            ch_versions = ch_versions.mix(DIFFCOEX.out.versions)
            ch_grn_results = ch_grn_results.mix(DIFFCOEX.out.grn_results)
        }

        if (methods.contains('grnboost2')) {
            GRNBOOST2(ch_metacells)
            ch_versions = ch_versions.mix(GRNBOOST2.out.versions)
            ch_grn_results = ch_grn_results.mix(GRNBOOST2.out.grn_results)
        }

        ch_grn_results.view()
    }

    
    emit:
    versions = ch_versions
}