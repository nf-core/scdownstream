include { SCANPY_PAGA            } from '../../modules/local/scanpy/paga'
include { SCANPY_RANKGENESGROUPS } from '../../modules/local/scanpy/rankgenesgroups'
include { LIANA_RANKAGGREGATE    } from '../../modules/local/liana/rankaggregate'

workflow PER_GROUP {
    take:
    ch_h5ad_with_neighbors // channel: [ meta, h5ad ], anndata objects with neighbors, one per embedding and annotation
    ch_h5ad_no_neighbors   // channel: [ meta, h5ad ], anndata objects without neighbors, one per annotation

    main:
    ch_versions      = channel.empty()
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()

    SCANPY_PAGA(ch_h5ad_with_neighbors)
    ch_versions      = ch_versions.mix(SCANPY_PAGA.out.versions)
    // ch_obsp       = ch_obsp.mix(SCANPY_PAGA.out.obsp)
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!params.skip_liana) {
        LIANA_RANKAGGREGATE(ch_h5ad_no_neighbors)
        ch_versions      = ch_versions.mix(LIANA_RANKAGGREGATE.out.versions)
        ch_uns           = ch_uns.mix(LIANA_RANKAGGREGATE.out.uns)
    }

    if (!params.skip_rankgenesgroups) {
        SCANPY_RANKGENESGROUPS(ch_h5ad_no_neighbors)
        ch_versions      = ch_versions.mix(SCANPY_RANKGENESGROUPS.out.versions)
        ch_uns           = ch_uns.mix(SCANPY_RANKGENESGROUPS.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(SCANPY_RANKGENESGROUPS.out.multiqc_files)
    }

    emit:
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions      // channel: [ versions.yml ]
}
