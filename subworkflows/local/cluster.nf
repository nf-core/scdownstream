include { SCANPY_NEIGHBORS } from '../../modules/local/scanpy/neighbors'
include { SCANPY_LEIDEN    } from '../../modules/local/scanpy/leiden'

workflow CLUSTER {
    take:
    ch_h5ad

    main:
    ch_versions = Channel.empty()
    ch_obs = Channel.empty()

    SCANPY_NEIGHBORS(ch_h5ad)
    ch_versions = ch_versions.mix(SCANPY_NEIGHBORS.out.versions)
    ch_h5ad = SCANPY_NEIGHBORS.out.h5ad

    ch_resolutions = Channel.from([0.5, 1.0])

    ch_h5ad = ch_h5ad.combine(ch_resolutions)
        .map{ meta, h5ad, resolution ->
                [meta + [resolution: resolution,
                        id: meta.integration + ":" + resolution],
            h5ad] }

    SCANPY_LEIDEN(ch_h5ad)
    ch_versions = ch_versions.mix(SCANPY_LEIDEN.out.versions)
    ch_obs = ch_obs.mix(SCANPY_LEIDEN.out.obs)

    emit:
    obs = ch_obs

    versions = ch_versions
}
