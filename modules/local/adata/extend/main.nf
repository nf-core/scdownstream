process ADATA_EXTEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45339bf761a2cf0cdb058492bc37f3df8b05b363731d491d1d3a14e9ba0b8f55/data'
        : 'community.wave.seqera.io/library/harmonypy_anndata_leidenalg_numpy_pruned:43066d5f86f18261'}"

    input:
    tuple (
        val(meta),
        path(base),
        path(obs, stageAs: 'obs/'),
        path(var, stageAs: 'var/'),
        path(obsm, stageAs: 'obsm/'),
        path(obsp, stageAs: 'obsp/'),
        path(uns, stageAs: 'uns/'),
        path(layers, stageAs: 'layers/')
    )

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv") , emit: metadata
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('extend.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_metadata.csv
    touch versions.yml
    """
}
