process DOUBLET_REMOVAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata_matplotlib_pandas_upsetplot:36786be0ea252945':
        'community.wave.seqera.io/library/anndata_matplotlib_pandas_upsetplot:f826c7c054eeef6d' }"

    input:
    tuple val(meta), path(h5ad), path(predictions)
    val(threshold)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*_mqc.json")             , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix    ?: "${meta.id}"
    template 'doublet_removal.py'

    stub:
    prefix = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
