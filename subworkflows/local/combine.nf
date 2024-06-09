include { ADATA_MERGE      } from '../../modules/local/adata/merge'
include { ADATA_UPSETGENES } from '../../modules/local/adata/upsetgenes'

workflow COMBINE {

    take:
    ch_h5ad // channel: [ val(meta), h5ad ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ADATA_MERGE(ch_h5ad.map { meta, h5ad -> h5ad }.collect()
                .map{ h5ads -> [[id: "merged"], h5ads] })
    ch_inner = ADATA_MERGE.out.inner
    ch_outer = ADATA_MERGE.out.outer
    ch_versions = ch_versions.mix(ADATA_MERGE.out.versions)

    ADATA_UPSETGENES(ch_outer)
    ch_versions = ch_versions.mix(ADATA_UPSETGENES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ADATA_UPSETGENES.out.multiqc_files)

    emit:

    multiqc_files = ch_multiqc_files
    versions      = ch_versions        // channel: [ versions.yml ]
}
