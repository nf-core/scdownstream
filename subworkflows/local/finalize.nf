include { ADATA_EXTEND } from '../../modules/local/adata/extend'
include { ADATA_TORDS  } from '../../modules/local/adata/tords'

workflow FINALIZE {
    take:
    ch_h5ad
    ch_obs
    ch_obsm

    main:
    ch_versions = Channel.empty()

    ADATA_EXTEND(ch_h5ad, ch_obs.flatten().collect(), ch_obsm.flatten().collect())
    ch_versions = ch_versions.mix(ADATA_EXTEND.out.versions)

    ADATA_TORDS(ADATA_EXTEND.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

    emit:
    versions = ch_versions
}
