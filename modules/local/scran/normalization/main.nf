process SCRAN_NORMALIZATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/671219a6661beb49512cf4183cf8e342435e85bfa37d70aa94fde6017c5b4988/data'
        : 'community.wave.seqera.io/library/bioconductor-scran_bioconductor-scuttle_bioconductor-scater_bioconductor-anndatar_pruned:d1de453471200309'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'normalization.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
