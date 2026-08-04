process ADATA_UNIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc259b0f0471f937af0396cdcb3834de97ca1fefa5563f160a217187600f69bc/data' :
        'community.wave.seqera.io/library/anndata_python_pyyaml:eec9fad329c5006d' }"

    input:
    tuple val(meta), path(h5ad)
    val batch_col
    val label_col
    val condition_col
    val unknown_label
    val symbol_col
    val counts_layer
    val duplicate_var_resolution
    val aggregate_isoforms
    val donor_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('unify.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
