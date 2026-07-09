include { CLUSTER_TARGETS                   } from '../cluster_targets'
include { SCANPY_NEIGHBORS as NEIGHBORS     } from '../../../modules/local/scanpy/neighbors'
include { SCANPY_LEIDEN as LEIDEN           } from '../../../modules/local/scanpy/leiden'
include { SCANPY_UMAP as UMAP               } from '../../../modules/local/scanpy/umap'
include { ADATA_ENTROPY as ENTROPY          } from '../../../modules/local/adata/entropy'
include { matchingAnalysisPlanRows          } from '../utils_nfcore_scdownstream_pipeline'
include { analysesFromPlanRows              } from '../utils_nfcore_scdownstream_pipeline'
include { deMethodsFromPlanRows             } from '../utils_nfcore_scdownstream_pipeline'

workflow CLUSTER {
    take:
    ch_input            // channel: [ meta, h5ad ]
    per_label           //   value: boolean
    global              //   value: boolean
    split_col           //   value: string
    analysis_plan_rows  //   value: list of plan rows; empty fields are wildcards
    default_resolutions //   value: list of resolution strings
    entropy_col         //   value: string
    embedding_key       //   value: string

    main:
    ch_obs = channel.empty()
    ch_obsm = channel.empty()
    ch_multiqc_files = channel.empty()

    CLUSTER_TARGETS (
        ch_input,
        per_label,
        global,
        split_col,
    )

    ch_h5ad = CLUSTER_TARGETS.out.targets.branch { meta, _h5ad ->
        has_neighbors: meta.integration == 'bbknn'
        needs_neighbors: true
    }

    NEIGHBORS (
        ch_h5ad.needs_neighbors,
        embedding_key
    )

    ch_h5ad_graph = NEIGHBORS.out.h5ad.mix(ch_h5ad.has_neighbors)

    UMAP (
        ch_h5ad_graph
    )
    ch_obsm = ch_obsm.mix(UMAP.out.obsm)

    ch_resolutions = channel.fromList(default_resolutions)

    ch_h5ad_for_leiden = UMAP.out.h5ad
        .combine(ch_resolutions)
        .map { meta, h5ad, resolution ->
            [matchingAnalysisPlanRows(analysis_plan_rows, meta, resolution), meta, h5ad, resolution]
        }
        .filter { matching_rows, _meta, _h5ad, _resolution ->
            !matching_rows.isEmpty()
        }
        .map { matching_rows, meta, h5ad, resolution ->
            [
                meta + [
                    resolution: resolution,
                    id: meta.id + '-' + resolution,
                ] + analysesFromPlanRows(matching_rows) + deMethodsFromPlanRows(matching_rows),
                h5ad,
            ]
        }

    ch_leiden = ch_h5ad_for_leiden.multiMap { meta, h5ad ->
        h5ad: [meta, h5ad]
        resolution: meta.resolution
        key_added: meta.id + '_leiden'
    }
    LEIDEN (
        ch_leiden.h5ad,
        ch_leiden.resolution,
        ch_leiden.key_added,
        true
    )
    ch_obs = ch_obs.mix(LEIDEN.out.obs)
    ch_h5ad_clustering = LEIDEN.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(LEIDEN.out.multiqc_files)

    ch_entropy = LEIDEN.out.h5ad
        .multiMap { meta, h5ad ->
            h5ad: [meta, h5ad]
            group_col: meta.id + '_leiden'
        }

    ENTROPY (
        ch_entropy.h5ad,
        ch_entropy.group_col,
        entropy_col
    )
    ch_obs = ch_obs.mix(ENTROPY.out.obs)
    ch_multiqc_files = ch_multiqc_files.mix(ENTROPY.out.multiqc_files)

    emit:
    obs             = ch_obs             // channel: [ pkl ]
    obsm            = ch_obsm            // channel: [ pkl ]
    h5ad_neighbors  = ch_h5ad_graph      // channel: [ meta, h5ad ]
    h5ad_clustering = ch_h5ad_clustering // channel: [ meta, h5ad ]
    multiqc_files   = ch_multiqc_files   // channel: [ json ]
}
