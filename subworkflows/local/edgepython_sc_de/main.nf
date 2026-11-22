include { EDGEPYTHON_SCDIFFERENTIAL } from '../../../modules/local/edgepython/sc_differential'
include { anndata                   } from 'plugin/nf-anndata'

workflow EDGEPYTHON_SC_DE {
    take:
    ch_h5ad // channel: [ meta, h5ad ] with meta.condition_col, meta.donor_col, and meta.obs_key
    reference_condition //   value: string

    main:
    ch_strata = ch_h5ad.flatMap { meta, h5ad ->
        def ad = anndata(h5ad)
        def celltype_col = meta.obs_key
        def condition_col = meta.condition_col
        def obs = ad.obs
        def celltypes = obs.get(celltype_col).data
        def conditions = obs.get(condition_col).data
        def donors = obs.get(meta.donor_col).data
        ad.close()

        def rowsByCelltype = [:].withDefault { [] }
        celltypes.eachWithIndex { celltype, index ->
            rowsByCelltype[celltype.toString()].add(index)
        }

        rowsByCelltype
            .findAll { celltype, rows ->
                def observedConditions = rows.collect { index -> conditions[index].toString() }.toSet()
                def observedDonors = rows.collect { index -> donors[index].toString() }.toSet()
                def eligible = rows.size() >= 10 && observedDonors.size() >= 2 && observedConditions.size() >= 2 && reference_condition in observedConditions
                if (!eligible) {
                    log.warn(
                        "Skipping edgepython_sc for '${meta.id}', ${celltype_col}=${celltype}. " + "Requires at least 10 cells, two donors and two conditions including '${reference_condition}'. " + "Found ${rows.size()} cells, ${observedDonors.size()} donors and conditions ${observedConditions.sort()}."
                    )
                }
                eligible
            }
            .collect { celltype, _rows ->
                [
                    meta + [celltype: celltype, condition_col: condition_col, de_method: 'edgepython_sc'],
                    h5ad,
                    celltype_col,
                    celltype as String,
                ]
            }
    }

    ch_edgepython = ch_strata.multiMap { meta, h5ad, celltype_col, celltype ->
        h5ad: [meta, h5ad]
        donor_col: meta.donor_col
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
