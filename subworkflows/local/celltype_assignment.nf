include { CELLTYPES_CELLTYPIST } from '../../modules/local/celltypes/celltypist'

workflow CELLTYPE_ASSIGNMENT {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()
    ch_obs = Channel.empty()

    if (params.celltypist_model) {
        celltypist_models = Channel.from(params.celltypist_model.split(','))

        CELLTYPES_CELLTYPIST(ch_h5ad, celltypist_models)
        ch_obs = ch_obs.mix(CELLTYPES_CELLTYPIST.out.obs)
        ch_versions = ch_versions.mix(CELLTYPES_CELLTYPIST.out.versions)
    }

    emit:
    obs = ch_obs

    versions = ch_versions
}
