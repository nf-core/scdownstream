process DECOUPLER_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5ca309507e02db5a619c02a9d8bb1d5b969ef6a815cceca6693784ee7fd061ec/data'
        : 'community.wave.seqera.io/library/decoupler_activity:8d8ed42929b0d74b'}"

    input:
    tuple val(meta), path(h5ad)
    val counts_layer
    val donor_col
    val celltype_col
    val condition_col
    val min_num_cells
    val min_total_counts

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}_samples.tsv", emit: samples
    path "*_mqc.json", emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('pseudobulk.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_samples.tsv
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
