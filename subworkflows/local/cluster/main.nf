include { ADATA_SPLITCOL as SPLITCOL    } from '../../../modules/local/adata/splitcol'
include { SCANPY_NEIGHBORS as NEIGHBORS } from '../../../modules/local/scanpy/neighbors'
include { SCANPY_LEIDEN as LEIDEN       } from '../../../modules/local/scanpy/leiden'
include { SCANPY_UMAP as UMAP           } from '../../../modules/local/scanpy/umap'
include { ADATA_ENTROPY as ENTROPY      } from '../../../modules/local/adata/entropy'
workflow CLUSTER {
    take:
    ch_input       // channel: [ integration, h5ad ]
    per_label      // value: boolean
    global         // value: boolean
    split_col      // value: string
    ch_resolutions // channel: [ string ]
    entropy_col    // value: string
    embedding_key  // value: string

    main:
    ch_versions = channel.empty()
    ch_obs = channel.empty()
    ch_obsm = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad = channel.empty()

    if (global) {
        ch_h5ad = ch_h5ad
            .mix(ch_input
                .map { meta, h5ad -> [meta + [subset: "global"], h5ad] })
    }

    if (per_label) {
        SPLITCOL (
            ch_input,
            split_col
        )
        ch_versions = ch_versions.mix(SPLITCOL.out.versions)

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
    ch_versions = ch_versions.mix(NEIGHBORS.out.versions)
    ch_h5ad = NEIGHBORS.out.h5ad.mix(ch_h5ad.has_neighbors)
    ch_h5ad_neighbours = NEIGHBORS.out.h5ad

    UMAP (
        ch_h5ad
    )
    ch_versions = ch_versions.mix(UMAP.out.versions)
    ch_obsm = ch_obsm.mix(UMAP.out.obsm)

    ch_h5ad = UMAP.out.h5ad
        .combine(ch_resolutions)
        .map { meta, h5ad, resolution ->
            [
                meta + [
                    resolution: resolution,
                    id: meta.integration + "-" + meta.subset + "-" + resolution,
                ],
                h5ad
            ]
        }

    ch_leiden = ch_h5ad.multiMap{ meta, h5ad ->
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
    ch_versions = ch_versions.mix(LEIDEN.out.versions)
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
    ch_versions = ch_versions.mix(ENTROPY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ENTROPY.out.multiqc_files)

    emit:
    obs             = ch_obs             // channel: [ pkl ]
    obsm            = ch_obsm            // channel: [ pkl ]
    h5ad_neighbors  = ch_h5ad_neighbours // channel: [ integration, h5ad ]
    h5ad_clustering = ch_h5ad_clustering // channel: [ integration, h5ad ]
    multiqc_files   = ch_multiqc_files   // channel: [ json ]
    versions        = ch_versions        // channel: [ versions.yml ]
}
