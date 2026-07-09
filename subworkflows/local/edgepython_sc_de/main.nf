include { EDGEPYTHON_SCDIFFERENTIAL } from '../../../modules/local/edgepython/sc_differential'
include { anndata       } from 'plugin/nf-anndata'

workflow EDGEPYTHON_SC_DE {
    take:
    ch_h5ad               // channel: [ meta, h5ad ] with meta.condition_col and meta.obs_key
    donor_col             //   value: string
    reference_condition   //   value: string

    main:
    ch_strata = ch_h5ad
        .map { meta, h5ad ->
            def ad = anndata(h5ad)
            def celltype_col = meta.obs_key
            def condition_col = meta.condition_col
            def celltypes = ad.obs[celltype_col].unique().toList()
            celltypes.collect { celltype ->
                [
                    meta + [celltype: celltype, condition_col: condition_col],
                    h5ad,
                    celltype_col,
                    celltype as String,
                ]
            }
        }
        .flatten()

    ch_edgepython = ch_strata.multiMap { meta, h5ad, celltype_col, celltype ->
        h5ad: [meta, h5ad]
        donor_col: donor_col
        condition_col: meta.condition_col
        celltype_col: celltype_col
        celltype_value: celltype
        reference_condition: reference_condition ?: ''
    }

    EDGEPYTHON_SCDIFFERENTIAL(
        ch_edgepython.h5ad,
        ch_edgepython.donor_col,
        ch_edgepython.condition_col,
        ch_edgepython.celltype_col,
        ch_edgepython.celltype_value,
        ch_edgepython.reference_condition,
    )

    emit:
    results = EDGEPYTHON_SCDIFFERENTIAL.out.results
}
