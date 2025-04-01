include { ADATA_PSEUDOBULKS } from '../../modules/local/pseudobulks'

workflow CREATE_PSEUDOBULKS {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()

    ADATA_PSEUDOBULKS(
            ch_h5ad.map{meta, h5ad -> [[id: 'pseudobulks'], h5ad]},
            params.pseudobulk_groups,
            params.pseudobulk_mode,
            params.pseudobulk_pseudoreplicates,
            params.pseudobulk_min_cells,
        )
    ch_versions = ch_versions.mix(ADATA_PSEUDOBULKS.out.versions)

    //ch_combined = ch_unfiltered.join(CELLBENDER_REMOVEBACKGROUND.out.barcodes)

    //ANNDATA_BARCODES(ch_combined)
    //ch_versions = ch_versions.mix(ANNDATA_BARCODES.out.versions)

    ch_h5ad = ADATA_PSEUDOBULKS.out.h5ad

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
