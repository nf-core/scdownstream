include { ADATA_EXTEND } from '../../modules/local/adata/extend'

workflow FINALIZE {
    take:
    ch_h5ad
    ch_obs
    ch_obsm

    main:
    ch_versions = Channel.empty()

    ADATA_EXTEND(ch_h5ad, ch_obs.flatten().collect(), ch_obsm.flatten().collect())
    ch_versions = ch_versions.mix(ADATA_EXTEND.out.versions)

    emit:
    versions = ch_versions
}
