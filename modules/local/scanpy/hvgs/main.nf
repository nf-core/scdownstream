process SCANPY_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/354081f1aea522edf9eb2503417bfd3ac227e5f6cc419192042d47a1457bfb15/data' :
        'community.wave.seqera.io/library/python_pyyaml_anndata_scanpy:7df6225e45ce62d7' }"

    input:
    tuple val(meta), path(h5ad)
    val n_hvgs
    path excluded_genes
    val input_layer
    val log_normalize
    val batch_key

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.pkl")                 , emit: var
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
    touch ${prefix}.pkl
    touch versions.yml
    """
}
