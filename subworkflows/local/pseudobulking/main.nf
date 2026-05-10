include { SCIMILARITY_PSEUDOBULK as PSEUDOBULK } from '../../../modules/local/scimilarity/pseudobulk'

workflow PSEUDOBULKING {
    take:
    ch_h5ad        // channel: [ integration, h5ad ]
    groupby_labels //   value: string
    min_num_cells  //   value: integer
    counts_layer   //   value: string

    main:
    ch_h5ad.view()
    ch_versions = channel.empty()

    PSEUDOBULK (
        ch_h5ad,
        counts_layer,
        groupby_labels,
        min_num_cells
    )
    ch_versions = ch_versions.mix(PSEUDOBULK.out.versions)
    ch_h5ad_pseudobulk = PSEUDOBULK.out.h5ad

    emit:
    h5ad_pseudobulk = ch_h5ad_pseudobulk // channel: [ integration, h5ad ]
    versions        = ch_versions        // channel: [ versions.yml ]
}
