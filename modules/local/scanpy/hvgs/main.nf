process SCANPY_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95e4cab1b5bbd9a3eb7a4139eac85d2440c6816ab8bdcede3944464ad5974843/data'
:         'community.wave.seqera.io/library/python_pyyaml_anndata_scanpy_pyarrow:803778b2933f5bc9' }"

    input:
    tuple val(meta), path(h5ad)
    val n_hvgs
    path excluded_genes
    val input_layer
    val log_normalize
    val batch_key

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.parquet")                 , emit: var
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('hvgs.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch versions.yml
    """
}
