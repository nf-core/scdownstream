include { SCANPY_RANKGENESGROUPS as RANKGENESGROUPS_CELLTYPES } from '../../modules/local/scanpy/rankgenesgroups'
include { SCANPY_RANKGENESGROUPS as RANKGENESGROUPS_CLUSTERS } from '../../modules/local/scanpy/rankgenesgroups'


workflow DIFFERENTIAL_EXPRESSION {
    take:
    ch_h5ad // channel: [ meta, h5ad ]
    ch_celltypes // channel: [ celltype ]
    ch_clusters // channel: [ meta, cluster ]

    main:
    ch_versions = Channel.empty()
    ch_uns = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_outdirs = Channel.empty()

    ch_input_celltypes = ch_h5ad.merge(ch_celltypes.flatten())
    RANKGENESGROUPS_CELLTYPES(ch_input_celltypes)
    ch_outdirs       = ch_outdirs.mix(RANKGENESGROUPS_CELLTYPES.out.outdir)
    ch_versions      = ch_versions.mix(RANKGENESGROUPS_CELLTYPES.out.versions)
    ch_uns           = ch_uns.mix(RANKGENESGROUPS_CELLTYPES.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(RANKGENESGROUPS_CELLTYPES.out.multiqc_files)

    ch_input_clusters = ch_h5ad
        .map { meta, h5ad -> [meta.integration, h5ad] }
        .combine(ch_clusters.map {meta, cluster -> [meta.integration, cluster]}, by: 0)
        .map { id, h5ad, cluster -> [[id: id], h5ad, cluster] }

    RANKGENESGROUPS_CLUSTERS(ch_input_clusters)
    ch_outdirs       = ch_outdirs.mix(RANKGENESGROUPS_CLUSTERS.out.outdir)
    ch_versions      = ch_versions.mix(RANKGENESGROUPS_CLUSTERS.out.versions)
    ch_uns           = ch_uns.mix(RANKGENESGROUPS_CLUSTERS.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(RANKGENESGROUPS_CLUSTERS.out.multiqc_files)

    emit:
    outdirs       = ch_outdirs       // channel: [ outdir ]
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions      // channel: [ versions.yml ]

}
