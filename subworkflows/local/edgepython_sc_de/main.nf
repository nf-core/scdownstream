include { EDGEPYTHON_SCDIFFERENTIAL } from '../../../modules/local/edgepython/sc_differential'
include { anndata                   } from 'plugin/nf-anndata'

workflow EDGEPYTHON_SC_DE {
    take:
    ch_h5ad // channel: [ meta, h5ad ] with meta.condition_col, meta.donor_col, and meta.obs_key
    reference_condition //   value: string

    main:
    ch_strata = ch_h5ad.flatMap { meta, h5ad ->
        def ad = anndata(h5ad)
        def obs = ad.obs
        def celltypes = obs.get(meta.obs_key).data
        def conditions = obs.get(meta.condition_col).data
        def donors = obs.get(meta.donor_col).data
        ad.close()

        def rowsByCelltype = (0..<celltypes.size()).groupBy { index -> celltypes[index].toString() }
        rowsByCelltype.findResults { celltype, rows ->
            def observedConditions = rows.collect { index -> conditions[index].toString() }.toSet()
            def donorCount = rows.collect { index -> donors[index].toString() }.toSet().size()
            if (rows.size() < 10 || donorCount < 2 || observedConditions.size() < 2 || !observedConditions.contains(reference_condition)) {
                log.warn("Skipping edgepython_sc for '${meta.id}', ${meta.obs_key}=${celltype}. Requires at least 10 cells, two donors and two conditions including '${reference_condition}'. Found ${rows.size()} cells, ${donorCount} donors and conditions ${observedConditions.sort()}.")
                return null
            }
            [meta + [celltype: celltype, de_method: 'edgepython_sc'], h5ad]
        }
    }

    ch_edgepython = ch_strata.multiMap { meta, h5ad ->
        h5ad: [meta, h5ad]
        donor_col: meta.donor_col
        condition_col: meta.condition_col
        celltype_col: meta.obs_key
        celltype_value: meta.celltype
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
