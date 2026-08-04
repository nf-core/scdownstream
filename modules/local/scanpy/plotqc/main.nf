process SCANPY_PLOTQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed34e104367da9d7fa40b39d78d86dce36948d95dbd2571abe923bc857bf192f/data' :
        'community.wave.seqera.io/library/python_pyyaml_scanpy:ec0559841ac88fb2' }"

    input:
    tuple val(meta), path(h5ad)
    val symbol_col
    path mito_genes
    val section_name
    val description

    output:
    tuple val(meta), path("*.png"), emit: plots
    path ("*_mqc.json")           , emit: multiqc_files
    path "versions.yml"           , emit: versions, topic: versions

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
