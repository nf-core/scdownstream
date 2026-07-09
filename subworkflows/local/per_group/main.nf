include { SCANPY_PAGA             } from '../../../modules/local/scanpy/paga'
include { LIANA_RANKAGGREGATE     } from '../../../modules/local/liana/rankaggregate'
include { DIFFERENTIAL_EXPRESSION } from '../differential_expression'
include { CYTETYPE                 } from '../../../modules/local/cytetype'
include { resolveDeMethodsWithPrerequisites } from '../utils_nfcore_scdownstream_pipeline'
include { cellLevelDeMethods                    } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkDeMethods                   } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsAnalysisEnabled        } from '../utils_nfcore_scdownstream_pipeline'

workflow PER_GROUP {
    take:
    ch_h5ad_with_neighbors        // channel: [ meta, h5ad ]
    ch_h5ad_no_neighbors          // channel: [ meta, h5ad ]
    skip_liana                    //   value: boolean
    cytetype_study_context        //   value: string
    de_methods_default            //   value: string
    pseudobulk_donor_col          //   value: string
    pseudobulk_min_num_cells      //   value: integer
    pseudobulk_min_total_counts   //   value: integer
    reference_condition           //   value: string

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad          = channel.empty()
    ch_obs           = channel.empty()

    SCANPY_PAGA(
        ch_h5ad_with_neighbors
            .filter { meta, _h5ad -> meta.analyses == null || 'paga' in meta.analyses }
    )
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!skip_liana) {
        LIANA_RANKAGGREGATE(
            ch_h5ad_no_neighbors
                .filter { meta, _h5ad -> meta.analyses == null || 'liana' in meta.analyses }
        )
        ch_uns = ch_uns.mix(LIANA_RANKAGGREGATE.out.uns)
    }

    ch_h5ad_for_de = ch_h5ad_no_neighbors
        .map { meta, h5ad ->
            [meta + resolveDeMethodsWithPrerequisites(meta, de_methods_default, cytetype_study_context), h5ad]
        }

    ch_h5ad_de = ch_h5ad_for_de
        .filter { meta, _h5ad ->
            def m = meta.de_methods_resolved
            m && (
                (rankGenesGroupsAnalysisEnabled(meta) && m.intersect(cellLevelDeMethods())) ||
                ((meta.analyses && 'pseudobulk_de' in meta.analyses) && m.intersect(pseudobulkDeMethods()))
            )
        }

    DIFFERENTIAL_EXPRESSION(
        ch_h5ad_de,
        pseudobulk_donor_col,
        pseudobulk_min_num_cells,
        pseudobulk_min_total_counts,
        reference_condition,
    )
    ch_uns           = ch_uns.mix(DIFFERENTIAL_EXPRESSION.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_EXPRESSION.out.multiqc_files)
    ch_h5ad          = ch_h5ad.mix(DIFFERENTIAL_EXPRESSION.out.h5ad)

    if (cytetype_study_context) {
        ch_cytetype = ch_h5ad
            .filter { meta, _h5ad -> meta.analyses == null || 'cytetype' in meta.analyses }
            .filter { meta, _h5ad ->
                meta.comparison_scope == 'global' && meta.de_method == 'wilcoxon'
            }
            .multiMap { meta, h5ad ->
                h5ad: [meta, h5ad]
                group_key: meta.obs_key
                rank_key: meta.rank_key
            }

        CYTETYPE(
            ch_cytetype.h5ad,
            "index",
            cytetype_study_context,
            ch_cytetype.group_key,
            ch_cytetype.rank_key,
        )
        ch_obs = ch_obs.mix(CYTETYPE.out.obs)
    }

    emit:
    uns           = ch_uns
    multiqc_files = ch_multiqc_files
    h5ad          = ch_h5ad
    obs           = ch_obs
}
