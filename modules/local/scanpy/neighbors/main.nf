process SCANPY_NEIGHBORS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c43892b5ea2991c98bd69e6da333b902a1e54635e24ca2cf3fd4c5972efb4384/data' :
        'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:e5a54bd2b6c9720d' }"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val(rep)
    val(n_pcs)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    template('neighbors.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
