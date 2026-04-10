process SCIMILARITY_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b07b44946319f0a77889ca315e4a48bef70c67bae06ce2603039a4995c75f0a/data'
        : 'community.wave.seqera.io/library/anndata_hnswlib_numcodecs_python_pruned:3f8ef15250e4fea7'}"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.pkl")                 , emit: obs
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    template('annotate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "versions.yml"
    """
}
