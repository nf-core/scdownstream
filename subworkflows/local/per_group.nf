include { SCANPY_PAGA            } from '../../modules/local/scanpy/paga'
include { SCANPY_RANKGENESGROUPS } from '../../modules/local/scanpy/rankgenesgroups'

workflow PER_GROUP {
    take:
    ch_h5ad

    main:
    ch_versions      = Channel.empty()
    ch_uns           = Channel.empty()
    ch_multiqc_files = Channel.empty()

    SCANPY_PAGA(ch_h5ad)
    ch_versions      = ch_versions.mix(SCANPY_PAGA.out.versions)
    // ch_obsp       = ch_obsp.mix(SCANPY_PAGA.out.obsp)
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!params.skip_rankgenesgroups) {
        SCANPY_RANKGENESGROUPS(ch_h5ad)
        ch_versions      = ch_versions.mix(SCANPY_RANKGENESGROUPS.out.versions)
        ch_uns           = ch_uns.mix(SCANPY_RANKGENESGROUPS.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(SCANPY_RANKGENESGROUPS.out.multiqc_files)
    }

    emit:
    uns           = ch_uns
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}