process SCANPY_UMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data' :
        'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf' }"

    input:
    tuple val(meta), path(h5ad, arity: 1)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl"                 , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('umap.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "X_${prefix}.pkl"
    touch "versions.yml"
    """
}
