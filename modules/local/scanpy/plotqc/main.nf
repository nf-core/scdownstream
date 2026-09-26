process SCANPY_PLOTQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf' }"

    input:
    tuple val(meta), path(h5ad)
    val symbol_col
    path mito_genes

    output:
    tuple val(meta), path("*.png"), emit: plots
    path ("*_mqc.json")           , emit: multiqc_files
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "QC Plots"
    description = task.ext.description ?: "Quality control plots"
    template('plotqc.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "QC Plots"
    description = task.ext.description ?: "Quality control plots"
    """
    touch ${prefix}_total_counts_vs_n_genes_by_counts.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
