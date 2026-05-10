process SCIMILARITY_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b07b44946319f0a77889ca315e4a48bef70c67bae06ce2603039a4995c75f0a/data'
        : 'community.wave.seqera.io/library/anndata_hnswlib_numcodecs_python_pruned:3f8ef15250e4fea7'}"

    input:
    tuple val(meta), path(h5ad)
    val(counts_layer)
    val(groupby_labels)
    val(min_num_cells)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    template('pseudobulk.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
