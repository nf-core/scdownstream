process SYMPHONY_HARMONYINTEGRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/512121548a21b4d1bb8acfd5e30a75c5c2103ddd00cf1de4713c682b7e6b5387/data'
            : 'community.wave.seqera.io/library/python_pyyaml_scanpy_pip_symphonypy:2198c27c5c9392d5'}"

    input:
    tuple val(meta), path(h5ad)
    val(batch_col)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl"                 , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    template('integrate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
