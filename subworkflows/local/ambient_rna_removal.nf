include { CELDA_DECONTX               } from '../../modules/local/celda/decontx'
include { CELLBENDER_REMOVEBACKGROUND } from '../../modules/local/cellbender/removebackground'
include { CELLBENDER_MERGE            } from '../../modules/local/cellbender/merge'
include { SOUPX                       } from '../../modules/local/soupx'
include { SCVITOOLS_SCAR              } from '../../modules/local/scvitools/scar'

workflow AMBIENT_RNA_REMOVAL {
    take:
    ch_pairing

    main:
    ch_versions = Channel.empty()

    if (params.ambient_removal == 'none') {
        println "AMBIENT_RNA_REMOVAL: Not performed since 'none' selected."
        ch_h5ad = ch_pairing.map{ meta, filtered, unfiltered -> [meta, filtered] }
    }
    else if (params.ambient_removal == 'decontx') {
        CELDA_DECONTX(ch_pairing)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = ch_versions.mix(CELDA_DECONTX.out.versions)
    }
    else if (params.ambient_removal == 'cellbender') {
        CELLBENDER_REMOVEBACKGROUND(ch_pairing.map{meta, filtered, unfiltered -> [meta, unfiltered]})
        ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

        CELLBENDER_MERGE(ch_pairing.map{ meta, filtered, raw -> [meta.id, meta, filtered, raw] }
            .join(CELLBENDER_REMOVEBACKGROUND.out.h5.map{ meta, h5 -> [meta.id, h5] }, by: 0, failOnMismatch: true)
            .map{ id, meta, filtered, raw, h5 -> [meta, filtered, raw, h5] })
        ch_h5ad = CELLBENDER_MERGE.out.h5ad
        ch_versions = ch_versions.mix(CELLBENDER_MERGE.out.versions)
    }
    else if (params.ambient_removal == 'soupx') {
        SOUPX(ch_pairing)
        ch_h5ad = SOUPX.out.h5ad
        ch_versions = ch_versions.mix(SOUPX.out.versions)
    }
    else if (params.ambient_removal == 'scar') {
        SCVITOOLS_SCAR(ch_pairing)
        ch_h5ad = SCVITOOLS_SCAR.out.h5ad
        ch_versions = SCVITOOLS_SCAR.out.versions
    }
    else {
        error "AMBIENT_RNA_REMOVAL: Unexpected value of param 'params.ambient_removal': '${params.ambient_removal}'."
    }

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
