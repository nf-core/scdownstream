process DECOUPLER_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63b778907daa8743d21b65174088943b8ea7cb6efaa79605c2bf2ba0b61d2a84/data'
        : 'community.wave.seqera.io/library/decoupler_pseudobulk:22ffb587b73301b8' }"

    input:
    tuple val(meta), path(h5ad)
    val(counts_layer)
    val(donor_col)
    val(celltype_col)
    val(condition_col)
    val(min_num_cells)
    val(min_total_counts)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}_samples.tsv"           , emit: samples
    path "versions.yml"                  , emit: versions, topic: versions

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
    touch versions.yml
    """
}
