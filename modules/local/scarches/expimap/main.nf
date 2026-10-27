process SCARCHES_EXPIMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6d/6d0046dded57ef25b360156cb77e1ddd34b95fae8a56b3d56171f6c12c3489c1/data'
        : 'community.wave.seqera.io/library/python_anndata_numpy_scipy_pruned:b66546c0fdf718b9'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    tuple val(meta2), path(reference_model)
    val(condition_col)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl", emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''
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
