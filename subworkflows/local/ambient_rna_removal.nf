include { CELDA_DECONTX         } from '../../modules/local/celda/decontx'
include { AMBIENTRNA_CELLBENDER } from '../../modules/local/ambientrna/cellbender'
include { ADATA_ADDH5           } from '../../modules/local/adata/addh5'

workflow AMBIENT_RNA_REMOVAL {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()

    if (params.ambient_removal == 'decontx') {
        CELDA_DECONTX(ch_h5ad)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = ch_versions.mix(CELDA_DECONTX.out.versions)
    }

    if (params.ambient_removal == 'cellbender') {
        AMBIENTRNA_CELLBENDER(ch_h5ad)
        ch_versions = ch_versions.mix(AMBIENTRNA_CELLBENDER.out.versions)

        ADATA_ADDH5(ch_h5ad.join(AMBIENTRNA_CELLBENDER.out.h5))
        ch_h5ad = ADATA_ADDH5.out.h5ad
        ch_versions = ch_versions.mix(ADATA_ADDH5.out.versions)
    }

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
