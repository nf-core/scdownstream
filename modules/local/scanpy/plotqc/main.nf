process SCANPY_PLOTQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c43892b5ea2991c98bd69e6da333b902a1e54635e24ca2cf3fd4c5972efb4384/data'
        : 'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:e5a54bd2b6c9720d'}"

    input:
    tuple val(meta), path(h5ad)
    val symbol_col
    path mito_genes
    val section_name

    output:
    tuple val(meta), path("*.png"), emit: plots
    path ("*_mqc.json"), emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('plotqc.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_total_counts_vs_n_genes_by_counts.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
