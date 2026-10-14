include { DECOUPLER_PSEUDOBULK as PSEUDOBULK } from '../../../modules/local/decoupler/pseudobulk'

workflow PSEUDOBULKING {
    take:
    ch_h5ad // channel: [ meta, anndata ] with meta.donor_col, meta.obs_key, meta.condition_col
    min_num_cells //   value: integer
    min_total_counts //   value: integer

    main:
    ch_pseudobulk = ch_h5ad.multiMap { meta, h5ad ->
        h5ad: [meta, h5ad]
        donor_col: meta.donor_col
        celltype_col: meta.obs_key
        condition_col: meta.condition_col
        counts_layer: meta.counts_layer ?: 'X'
        min_num_cells: min_num_cells
        min_total_counts: min_total_counts
    }

    PSEUDOBULK(
        ch_pseudobulk.h5ad,
        ch_pseudobulk.counts_layer,
        ch_pseudobulk.donor_col,
        ch_pseudobulk.celltype_col,
        ch_pseudobulk.condition_col,
        ch_pseudobulk.min_num_cells,
        ch_pseudobulk.min_total_counts,
        'zarr',
    )

    emit:
    anndata = PSEUDOBULK.out.anndata
    samples = PSEUDOBULK.out.samples
    multiqc_files = PSEUDOBULK.out.multiqc_files
}
