include { SCANPY_PAGA                       } from '../../../modules/local/scanpy/paga'
include { LIANA_RANKAGGREGATE               } from '../../../modules/local/liana/rankaggregate'
include { DIFFERENTIAL_EXPRESSION           } from '../differential_expression'
include { CLUSTER_ANNOTATION                } from '../cluster_annotation'
include { resolveDeMethodsWithPrerequisites } from '../utils_nfcore_scdownstream_pipeline'
include { cellLevelDeMethods                } from '../utils_nfcore_scdownstream_pipeline'
include { pseudobulkingEnabled              } from '../utils_nfcore_scdownstream_pipeline'
include { rankGenesGroupsAnalysisEnabled    } from '../utils_nfcore_scdownstream_pipeline'

workflow PER_GROUP {
    take:
    ch_h5ad_with_neighbors         // channel: [ meta, h5ad ]
    ch_h5ad_no_neighbors           // channel: [ meta, h5ad ]
    skip_liana                     //   value: boolean
    liana_n_perms                  //   value: integer
    liana_max_cells                //   value: integer or null
    liana_subsample_strategy       //   value: string
    liana_subsample_seed           //   value: integer
    cytetype_study_context         //   value: string
    de_methods_default             //   value: string
    pseudobulk                     //   value: boolean
    pseudobulk_min_num_cells       //   value: integer
    pseudobulk_min_total_counts    //   value: integer
    ch_per_cell_annotation_columns // channel: string
    reference_condition            //   value: string
    interesting_genes              //   value: string (path) or []

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad          = channel.empty()

    SCANPY_PAGA(
        ch_h5ad_with_neighbors
            .filter { meta, _h5ad -> meta.analyses == null || 'paga' in meta.analyses }
    )
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!skip_liana) {
        LIANA_RANKAGGREGATE(
            ch_h5ad_no_neighbors
                .filter { meta, _h5ad -> meta.analyses == null || 'liana' in meta.analyses },
            liana_n_perms,
            liana_max_cells ?: 0,
            liana_subsample_strategy,
            liana_subsample_seed,
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
                pseudobulkingEnabled(meta, pseudobulk) ||
                ((meta.analyses == null || 'de' in meta.analyses) && 'edgepython_sc' in m)
            )
        }

    DIFFERENTIAL_EXPRESSION(
        ch_h5ad_de,
        pseudobulk,
        pseudobulk_min_num_cells,
        pseudobulk_min_total_counts,
        reference_condition,
        interesting_genes,
    )
    ch_uns           = ch_uns.mix(DIFFERENTIAL_EXPRESSION.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_EXPRESSION.out.multiqc_files)
    ch_h5ad          = ch_h5ad.mix(DIFFERENTIAL_EXPRESSION.out.h5ad)

    CLUSTER_ANNOTATION(
        ch_h5ad_no_neighbors,
        ch_h5ad,
        ch_per_cell_annotation_columns,
        cytetype_study_context,
    )
    ch_obs = CLUSTER_ANNOTATION.out.obs

    emit:
    uns           = ch_uns
    multiqc_files = ch_multiqc_files.flatten()
    h5ad          = ch_h5ad
    obs           = ch_obs
}
