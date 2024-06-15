include { CELDA_DECONTX         } from '../../modules/local/celda/decontx'
include { AMBIENTRNA_CELLBENDER } from '../../modules/local/ambientrna/cellbender'

workflow AMBIENT_RNA_REMOVAL {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()

    if (params.ambient_removal == 'decontx') {
        CELDA_DECONTX(ch_h5ad)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = CELDA_DECONTX.out.versions
    }

    if (params.ambient_removal == 'cellbender') {
        AMBIENTRNA_CELLBENDER(ch_h5ad)
        ch_h5ad = AMBIENTRNA_CELLBENDER.out.h5ad
        ch_versions = AMBIENTRNA_CELLBENDER.out.versions
    }

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
