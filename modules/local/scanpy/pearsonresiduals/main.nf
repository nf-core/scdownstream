process SCANPY_PEARSONRESIDUALS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f692c0c793f5ac7da42963b22084153c26be689e512864095832318e77dee2f9/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy:641b33efb395c24e' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "pearson_residuals.np*"           , emit: layers
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('pearsonresiduals.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch pearson_residuals.npz
    touch versions.yml
    """
}
