process SCDS {
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/nicotru/scds:7788dbeb87bc7eec"

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
    template('scds.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv
    touch versions.yml
    """
}
