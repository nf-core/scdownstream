include { CUSTOM_AGGREGATE_PER_CELL_ANNOTATIONS } from '../../../modules/local/custom/aggregate_per_cell_annotations'
include { CYTETYPE                          } from '../../../modules/local/cytetype'
include { analysisEnabled                   } from '../utils_nfcore_scdownstream_pipeline'

workflow CLUSTER_ANNOTATION {
    take:
    ch_h5ad                        // channel: [ meta, h5ad ]
    ch_h5ad_de                     // channel: [ meta, h5ad ] — global DE results for CyteType
    ch_per_cell_annotation_columns // channel: string
    cytetype_study_context         //   value: string

    main:
    ch_obs = channel.empty()

    ch_aggregate_per_cell_annotations = ch_h5ad
        .filter { meta, _h5ad -> analysisEnabled(meta, 'aggregate_per_cell_annotation') }
        .combine(ch_per_cell_annotation_columns)
        .map { meta, h5ad, annotation_col ->
            [meta + [id: "${meta.id}:${annotation_col}"], h5ad, meta.obs_key, annotation_col, meta.integration ?: 'unknown', meta.resolution ?: 'unknown']
        }
        .multiMap { meta, h5ad, cluster_col, annotation_col, integration, resolution ->
            h5ad: [meta, h5ad]
            cluster_col: cluster_col
            annotation_col: annotation_col
            integration: integration
            resolution: resolution
        }

    CUSTOM_AGGREGATE_PER_CELL_ANNOTATIONS(
        ch_aggregate_per_cell_annotations.h5ad,
        ch_aggregate_per_cell_annotations.cluster_col,
        ch_aggregate_per_cell_annotations.annotation_col,
        ch_aggregate_per_cell_annotations.integration,
        ch_aggregate_per_cell_annotations.resolution,
    )
    ch_obs = ch_obs.mix(CUSTOM_AGGREGATE_PER_CELL_ANNOTATIONS.out.obs)

    if (cytetype_study_context) {
        ch_cytetype = ch_h5ad_de
            .filter { meta, _h5ad -> analysisEnabled(meta, 'cytetype') }
            .filter { meta, _h5ad ->
                meta.comparison_scope == 'global' && meta.de_method == 'wilcoxon'
            }
            .multiMap { meta, h5ad ->
                h5ad: [meta, h5ad]
                group_key: meta.obs_key
                rank_key: meta.rank_key
                integration: meta.integration ?: 'unknown'
                resolution: meta.resolution ?: 'unknown'
            }

        CYTETYPE(
            ch_cytetype.h5ad,
            "index",
            cytetype_study_context,
            ch_cytetype.group_key,
            ch_cytetype.rank_key,
            ch_cytetype.integration,
            ch_cytetype.resolution,
        )
        ch_obs = ch_obs.mix(CYTETYPE.out.obs)
    }

    emit:
    obs = ch_obs
}
