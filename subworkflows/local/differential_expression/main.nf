include { PSEUDOBULKING       } from '../pseudobulking'
include { PSEUDOBULK_DE        } from '../pseudobulk_de'
include { RANK_GENES_GROUPS    } from '../rank_genes_groups'
include { EDGEPYTHON_SC_DE     } from '../edgepython_sc_de'
include { rankGenesGroupsMethods         } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsAnalysisEnabled } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeMethods            } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkingEnabled           } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeEnabled            } from '../utils_nfcore_scdownstream_pipeline'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    ch_h5ad                     // channel: [ meta, h5ad ] with de_methods_resolved, obs_key, condition_col, analyses
    pseudobulk                  //   value: boolean
    pseudobulk_min_num_cells    //   value: integer
    pseudobulk_min_total_counts //   value: integer
    reference_condition         //   value: string
    interesting_genes           //   value: string (path) or []

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad_out      = channel.empty()

    ch_h5ad_pseudobulk = ch_h5ad
        .filter { meta, _h5ad -> pseudobulkingEnabled(meta, pseudobulk) }

    PSEUDOBULKING(
        ch_h5ad_pseudobulk,
        pseudobulk_min_num_cells,
        pseudobulk_min_total_counts,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PSEUDOBULKING.out.multiqc_files)

    PSEUDOBULK_DE(
        PSEUDOBULKING.out.h5ad
            .filter { meta, _h5ad -> pseudobulkDeEnabled(meta) },
        reference_condition,
        interesting_genes,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PSEUDOBULK_DE.out.multiqc_files)

    ch_h5ad_rgg = ch_h5ad
        .filter { meta, _h5ad ->
            rankGenesGroupsAnalysisEnabled(meta) &&
                meta.de_methods_resolved.intersect(rankGenesGroupsMethods())
        }

    RANK_GENES_GROUPS(ch_h5ad_rgg, interesting_genes)
    ch_uns           = ch_uns.mix(RANK_GENES_GROUPS.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(RANK_GENES_GROUPS.out.multiqc_files)
    ch_h5ad_out      = ch_h5ad_out.mix(RANK_GENES_GROUPS.out.h5ad)

    EDGEPYTHON_SC_DE(
        ch_h5ad
            .filter { meta, _h5ad ->
                (meta.analyses == null || 'de' in meta.analyses) &&
                    'edgepython_sc' in meta.de_methods_resolved
            },
        reference_condition,
        interesting_genes,
    )
    ch_multiqc_files = ch_multiqc_files.mix(EDGEPYTHON_SC_DE.out.multiqc_files)

    emit:
    uns           = ch_uns
    multiqc_files = ch_multiqc_files
    h5ad          = ch_h5ad_out
}
