include { SCANPY_PAGA            } from '../../modules/local/scanpy/paga'
include { SCANPY_RANKGENESGROUPS } from '../../modules/local/scanpy/rankgenesgroups'
include { LIANA_RANKAGGREGATE    } from '../../modules/local/liana/rankaggregate'

workflow PER_GROUP {
    take:
    ch_h5ad_both           // channel: [ integration, h5ad ]
    ch_h5ad_with_neighbors // channel: [ integration, h5ad ]
    ch_h5ad_no_neighbors   // channel: [ merged, h5ad ]

    main:
    ch_versions      = Channel.empty()
    ch_uns           = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_with_neighbors = ch_h5ad_both.mix(ch_h5ad_with_neighbors)
    ch_no_neighbors   = ch_h5ad_both.mix(ch_h5ad_no_neighbors)

    SCANPY_PAGA(ch_with_neighbors)
    ch_versions      = ch_versions.mix(SCANPY_PAGA.out.versions)
    // ch_obsp       = ch_obsp.mix(SCANPY_PAGA.out.obsp)
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!params.skip_liana) {
        LIANA_RANKAGGREGATE(ch_no_neighbors)
        ch_versions      = ch_versions.mix(LIANA_RANKAGGREGATE.out.versions)
        ch_uns           = ch_uns.mix(LIANA_RANKAGGREGATE.out.uns)
    }

    emit:
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions      // channel: [ versions.yml ]
}