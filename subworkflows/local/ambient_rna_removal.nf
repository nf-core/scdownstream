include { CELDA_DECONTX } from '../../modules/local/celda/decontx'

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

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
