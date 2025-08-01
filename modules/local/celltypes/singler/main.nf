process CELLTYPES_SINGLER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/saditya88/singler:0.0.1'


    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    tuple val(meta2), val(names), val(labels), path(references)

    output:
    //tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv")             , emit: obs
    tuple val(meta), path("*_distribution.pdf"), emit: distribution
    tuple val(meta), path("*_heatmap.pdf")     , emit: heatmap
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SINGLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    template 'singleR.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SINGLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    touch ${prefix}_distribution.pdf
    touch ${prefix}_heatmap.pdf
    touch ${prefix}.csv
    touch versions.yml
    """
}
