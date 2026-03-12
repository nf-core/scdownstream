process SCDBLFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-biocparallel_bioconductor-rhdf5_bioconductor-scdblfinder_pruned:28444625ead66428' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-biocparallel_bioconductor-rhdf5_bioconductor-scdblfinder_pruned:0f9db6b0855861de' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.csv"), emit: predictions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scdblfinder.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv
    touch versions.yml
    """
}
