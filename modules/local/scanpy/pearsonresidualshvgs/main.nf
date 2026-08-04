process SCANPY_PEARSONRESIDUALS_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed34e104367da9d7fa40b39d78d86dce36948d95dbd2571abe923bc857bf192f/data' :
        'community.wave.seqera.io/library/python_pyyaml_scanpy:ec0559841ac88fb2' }"

    input:
    tuple val(meta), path(h5ad)
    val n_genes
    path excluded_genes
    val(batch_key)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.pkl")                 , emit: var
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('pearsonresidualshvgs.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch versions.yml
    """
}
