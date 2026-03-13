process SCARCHES_EXPIMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
            ? 'https://wave.seqera.io/view/builds/bd-e532922f69f9a648_1'
            : 'community.wave.seqera.io/library/pip_scarches:7e8c7e577326f6ad'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    tuple val(meta2), path(reference_model)
    val(batch_col)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl", emit: obsm
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    template('expimap.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
