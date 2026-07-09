include { PSEUDOBULKING       } from '../pseudobulking'
include { PSEUDOBULK_DE        } from '../pseudobulk_de'
include { RANK_GENES_GROUPS    } from '../rank_genes_groups'
include { EDGEPYTHON_SC_DE     } from '../edgepython_sc_de'
include { rankGenesGroupsMethods         } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsAnalysisEnabled } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeMethods    } from '../utils_nfcore_scdownstream_pipeline'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    ch_h5ad // channel: [ meta, h5ad ] with de_methods_resolved, obs_key, condition_col, analyses
    pseudobulk_donor_col        //   value: string
    pseudobulk_min_num_cells    //   value: integer
    pseudobulk_min_total_counts //   value: integer
    reference_condition         //   value: string

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad_out      = channel.empty()

    ch_h5ad_pseudobulk_de = ch_h5ad
        .filter { meta, _h5ad ->
            meta.analyses && 'pseudobulk_de' in meta.analyses &&
                meta.de_methods_resolved.intersect(pseudobulkDeMethods())
        }

    PSEUDOBULKING(
        ch_h5ad_pseudobulk_de,
        pseudobulk_donor_col,
        pseudobulk_min_num_cells,
        pseudobulk_min_total_counts,
    )

    PSEUDOBULK_DE(
        PSEUDOBULKING.out.h5ad,
        reference_condition,
    )

    ch_h5ad_rgg = ch_h5ad
        .filter { meta, _h5ad ->
            rankGenesGroupsAnalysisEnabled(meta) &&
                meta.de_methods_resolved.intersect(rankGenesGroupsMethods())
        }

    RANK_GENES_GROUPS(ch_h5ad_rgg)
    ch_uns           = ch_uns.mix(RANK_GENES_GROUPS.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(RANK_GENES_GROUPS.out.multiqc_files)
    ch_h5ad_out      = ch_h5ad_out.mix(RANK_GENES_GROUPS.out.h5ad)

    EDGEPYTHON_SC_DE(
        ch_h5ad
            .filter { meta, _h5ad ->
                (meta.analyses == null || 'de' in meta.analyses) &&
                    'edgepython_sc' in meta.de_methods_resolved
            },
        pseudobulk_donor_col,
        reference_condition,
    )

    emit:
    uns           = ch_uns
    multiqc_files = ch_multiqc_files
    h5ad          = ch_h5ad_out
}
