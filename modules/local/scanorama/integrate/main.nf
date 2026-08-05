process SCANORAMA_INTEGRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/345eb2fb37907b25b7b18f33c7a933799874b9177a89abdff6478b28f414cdcc/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_scanorama_pyarrow:100d7fce0082a614' }"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val(batch_col)
    val(input_layer)
    val(log_normalize)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.parquet"                 , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('integrate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "X_${prefix}.parquet"
    touch "versions.yml"
    """
}
