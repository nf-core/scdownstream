include { SCANPY_NEIGHBORS } from '../../modules/local/scanpy/neighbors'
include { SCANPY_LEIDEN    } from '../../modules/local/scanpy/leiden'
include { SCANPY_UMAP      } from '../../modules/local/scanpy/umap'
include { SCANPY_PAGA      } from '../../modules/local/scanpy/paga'

workflow CLUSTER {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()
    ch_obs = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obsp = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_h5ad.branch{ meta, h5ad ->
        has_neighbors: meta.integration == "bbknn"
        needs_neighbors: true
    }.set { ch_h5ad }

    SCANPY_NEIGHBORS(ch_h5ad.needs_neighbors)
    ch_versions = ch_versions.mix(SCANPY_NEIGHBORS.out.versions)
    ch_h5ad = SCANPY_NEIGHBORS.out.h5ad.mix(ch_h5ad.has_neighbors)

    SCANPY_UMAP(ch_h5ad)
    ch_versions = ch_versions.mix(SCANPY_UMAP.out.versions)
    ch_obsm = ch_obsm.mix(SCANPY_UMAP.out.obsm)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_UMAP.out.multiqc_files)

    ch_resolutions = Channel.from(params.clustering_resolutions.split(","))

    ch_h5ad = SCANPY_UMAP.out.h5ad.combine(ch_resolutions)
        .map{ meta, h5ad, resolution ->
                [meta + [resolution: resolution,
                        id: meta.integration + "-" + resolution],
            h5ad] }

    SCANPY_LEIDEN(ch_h5ad)
    ch_versions = ch_versions.mix(SCANPY_LEIDEN.out.versions)
    ch_obs = ch_obs.mix(SCANPY_LEIDEN.out.obs)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_LEIDEN.out.multiqc_files)

    SCANPY_PAGA(SCANPY_LEIDEN.out.h5ad)
    ch_versions = ch_versions.mix(SCANPY_PAGA.out.versions)
    ch_obsp = ch_obsp.mix(SCANPY_PAGA.out.obsp)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    emit:
    obs = ch_obs
    obsm = ch_obsm
    obsp = ch_obsp

    multiqc_files = ch_multiqc_files
    versions = ch_versions
}
