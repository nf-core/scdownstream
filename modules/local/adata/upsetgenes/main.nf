process ADATA_UPSETGENES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/373aa42c5a49c2b6892606e91b80eb528c90d588261e51b1272071d0e608c621/data' :
        'community.wave.seqera.io/library/adata_upsetgenes:99308af0235c51d3' }"

    input:
    tuple val(meta), val(names), path(h5ads)

    output:
    tuple val(meta), path("*.png"), emit: plot, optional: true
    path ("*_mqc.json")           , emit: multiqc_files, optional: true
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('upsetplot.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
