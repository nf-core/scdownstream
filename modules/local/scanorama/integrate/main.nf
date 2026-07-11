process SCANORAMA_INTEGRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/ee5ffb0917263d09101dd05dc45c76aa239af97c58f2f5a368ac8a67a3eede3c/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_scanorama:095181a9ed452380' }"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val(batch_col)
    val(input_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl"                 , emit: obsm
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
    touch "X_${prefix}.pkl"
    touch "versions.yml"
    """
}
