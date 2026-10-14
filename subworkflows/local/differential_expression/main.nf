include { PSEUDOBULKING } from '../pseudobulking'
include { PSEUDOBULK_DE } from '../pseudobulk_de'
include { RANK_GENES_GROUPS } from '../rank_genes_groups'
include { EDGEPYTHON_SC_DE } from '../edgepython_sc_de'
include { CUSTOM_VOLCANOPLOT } from '../../../modules/local/custom/volcanoplot'
include { rankGenesGroupsMethods } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsAnalysisEnabled } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeMethods } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkingEnabled } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeEnabled } from '../utils_nfcore_scdownstream_pipeline'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    ch_h5ad // channel: [ meta, anndata ] with de_methods_resolved, obs_key, condition_col, analyses
    pseudobulk //   value: boolean
    decoupler //   value: boolean
    pseudobulk_min_num_cells //   value: integer
    pseudobulk_min_total_counts //   value: integer
    reference_condition //   value: string
    interesting_genes //   value: string (path) or []
    save_intermediates //   value: boolean
    cytetype_study_context //   value: string

    main:
    ch_uns = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad_out = channel.empty()

    ch_h5ad_pseudobulk = ch_h5ad.filter { meta, _h5ad -> pseudobulkingEnabled(meta, pseudobulk, decoupler) }

    PSEUDOBULKING(
        ch_h5ad_pseudobulk,
        pseudobulk_min_num_cells,
        pseudobulk_min_total_counts,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PSEUDOBULKING.out.multiqc_files)

    PSEUDOBULK_DE(
        PSEUDOBULKING.out.anndata.filter { meta, _h5ad -> pseudobulkDeEnabled(meta) },
        reference_condition,
    )

    ch_h5ad_rgg = ch_h5ad.filter { meta, _h5ad ->
        rankGenesGroupsAnalysisEnabled(meta) && meta.de_methods_resolved.intersect(rankGenesGroupsMethods())
    }

    RANK_GENES_GROUPS(ch_h5ad_rgg, save_intermediates, cytetype_study_context)
    ch_uns = ch_uns.mix(RANK_GENES_GROUPS.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(RANK_GENES_GROUPS.out.multiqc_files)
    ch_h5ad_out = ch_h5ad_out.mix(RANK_GENES_GROUPS.out.anndata)

    EDGEPYTHON_SC_DE(
        ch_h5ad.filter { meta, _h5ad ->
            (meta.analyses == null || 'de' in meta.analyses) && 'edgepython_sc' in meta.de_methods_resolved
        },
        reference_condition,
    )

    ch_de_results = RANK_GENES_GROUPS.out.results
        .mix(PSEUDOBULK_DE.out.results)
        .mix(EDGEPYTHON_SC_DE.out.results)
        .flatMap { meta, files ->
            def collected = files instanceof Collection ? files : [files]
            collected.collect { pq -> [meta, pq] }
        }

    CUSTOM_VOLCANOPLOT(ch_de_results, interesting_genes ?: [])
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_VOLCANOPLOT.out.multiqc_files)

    emit:
    uns = ch_uns
    multiqc_files = ch_multiqc_files
    anndata = ch_h5ad_out
    pseudobulk_anndata = PSEUDOBULKING.out.anndata
    pseudobulk_de_results = PSEUDOBULK_DE.out.results
}
