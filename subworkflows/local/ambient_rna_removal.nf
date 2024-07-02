include { CELDA_DECONTX               } from '../../modules/local/celda/decontx'
include { CELLBENDER_REMOVEBACKGROUND } from '../../modules/local/cellbender/removebackground'
include { CELLBENDER_MERGE            } from '../../modules/local/cellbender/merge'
include { SOUPX                       } from '../../modules/local/soupx'
include { SCVITOOLS_SCAR              } from '../../modules/local/scvitools/scar'

workflow AMBIENT_RNA_REMOVAL {
    take:
    ch_filtered
    ch_raw

    main:
    ch_versions = Channel.empty()

    if (params.ambient_removal == 'none') {
        println "AMBIENT_RNA_REMOVAL: Not performed since 'none' selected."
    }
    else if (params.ambient_removal == 'decontx') {
        CELDA_DECONTX(ch_filtered)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = ch_versions.mix(CELDA_DECONTX.out.versions)
    }
    else if (params.ambient_removal == 'cellbender') {
        CELLBENDER_REMOVEBACKGROUND(ch_raw)
        ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

        CELLBENDER_MERGE(ch_filtered.map{ meta, h5ad -> [meta.id, meta, h5ad] }
            .join(ch_raw.join(CELLBENDER_REMOVEBACKGROUND.out.h5)
                .map{ meta, raw, cellbender -> [meta.id, raw, cellbender] }, failOnMismatch: true)
            .map{ id, meta, filtered, raw, cellbender -> [meta, filtered, raw, cellbender] }
            )
        ch_h5ad = CELLBENDER_MERGE.out.h5ad
        ch_versions = ch_versions.mix(CELLBENDER_MERGE.out.versions)
    }
    else if (params.ambient_removal == 'soupx') {
        SOUPX(ch_filtered.map{ meta, h5ad -> [meta.id, meta, h5ad]}
            .join(ch_raw.map{ meta, h5ad -> [meta.id, h5ad]}, failOnMismatch: true)
            .map{ id, meta, h5ad, raw -> [meta, h5ad, raw] })
        ch_h5ad = SOUPX.out.h5ad
        ch_versions = ch_versions.mix(SOUPX.out.versions)
    }
    else if (params.ambient_removal == 'scar') {
        SCVITOOLS_SCAR(ch_filtered.map{ meta, h5ad -> [meta.id, meta, h5ad]}
            .join(ch_raw.map{ meta, h5ad -> [meta.id, h5ad]}, failOnMismatch: true)
            .map{ id, meta, h5ad, raw -> [meta, h5ad, raw] })
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
