process SYMPHONY_MAPEMBEDDING {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fef1581af3ad5950a2e168b8bf07dacefd670a2877a53d0b16196a4ce3f7ac56/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_pip_pruned:b85ec8f003db501a' }"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(reference_h5ad, stageAs: 'reference/reference.h5ad')
    val(batch_col)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.parquet"                 , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    template('map_embedding.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
