process SCANPY_NEIGHBORS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f36c00fcbaa1da1aaa6cf25cd8a562efb6728457d4a54b46fa06395aebd2083e/data' :
        'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:53120062012176bc' }"

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
