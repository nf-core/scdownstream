include { ADATA_SPLITCOL as SPLITCOL    } from '../../../modules/local/adata/splitcol'
include { SCANPY_NEIGHBORS as NEIGHBORS } from '../../../modules/local/scanpy/neighbors'
include { SCANPY_LEIDEN as LEIDEN       } from '../../../modules/local/scanpy/leiden'
include { SCANPY_UMAP as UMAP           } from '../../../modules/local/scanpy/umap'
include { ADATA_ENTROPY as ENTROPY      } from '../../../modules/local/adata/entropy'
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
    ch_h5ad = channel.empty()

    ch_input_by_subset = ch_input.branch { meta, _h5ad ->
        already_split: meta.subset != null
        needs_split: true
    }

    ch_h5ad = ch_h5ad.mix(
        ch_input_by_subset.already_split
            .map { meta, h5ad -> [meta + [already_split: true], h5ad] }
    )

    if (global) {
        ch_h5ad = ch_h5ad
            .mix(ch_input_by_subset.needs_split
                .map { meta, h5ad -> [meta + [subset: "global"], h5ad] })
    }

    if (per_label) {
        SPLITCOL (
            ch_input_by_subset.needs_split,
            split_col
        )

        ch_h5ad = ch_h5ad.mix(
            SPLITCOL.out.h5ad
                .transpose()
                .map { meta, h5ad -> [meta + [subset: h5ad.simpleName], h5ad] }
        )
    }

    ch_h5ad = ch_h5ad
        .map {
            meta, h5ad ->
            [meta + [id: meta.integration + "-" + meta.subset], h5ad]
        }

    ch_h5ad = ch_h5ad.branch { meta, _h5ad ->
        has_neighbors: meta.integration == "bbknn"
        needs_neighbors: true
    }

    NEIGHBORS (
        ch_h5ad.needs_neighbors,
        embedding_key
    )
    ch_h5ad = NEIGHBORS.out.h5ad.mix(ch_h5ad.has_neighbors)
    ch_h5ad_neighbours = NEIGHBORS.out.h5ad

    ch_h5ad_for_umap = ch_h5ad
        .map { meta, h5ad ->
            meta.already_split
                ? [meta + [id: meta.id + "-umap", cluster_id: meta.id], h5ad]
                : [meta, h5ad]
        }

    UMAP (
        ch_h5ad_for_umap
    )
    ch_obsm = ch_obsm.mix(UMAP.out.obsm)

    ch_resolutions = channel.fromList(default_resolutions)

    ch_h5ad_for_leiden = UMAP.out.h5ad
        .map { meta, h5ad ->
            meta.cluster_id
                ? [meta.findAll { key, _value -> key != 'cluster_id' } + [id: meta.cluster_id], h5ad]
                : [meta, h5ad]
        }
        .combine(ch_resolutions)
        .filter { meta, _h5ad, resolution ->
            analysis_plan_rows.any { row ->
                (!row.integration || row.integration == meta.integration) &&
                (!row.subset || row.subset == meta.subset) &&
                (!row.resolution || (row.resolution as String) == resolution)
            }
        }
        .map { meta, h5ad, resolution ->
            [
                meta + [
                    resolution: resolution,
                    id: meta.integration + "-" + meta.subset + "-" + resolution,
                ],
                h5ad,
            ]
        }

    ch_leiden = ch_h5ad_for_leiden.multiMap{ meta, h5ad ->
        h5ad: [meta, h5ad]
        resolution: meta.resolution
        key_added: meta.id + "_leiden"
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
            group_col: meta.id + "_leiden"
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
    h5ad_neighbors  = ch_h5ad_neighbours // channel: [ integration, h5ad ]
    h5ad_clustering = ch_h5ad_clustering // channel: [ integration, h5ad ]
    multiqc_files   = ch_multiqc_files   // channel: [ json ]
}
