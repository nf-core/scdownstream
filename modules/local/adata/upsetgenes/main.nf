process ADATA_UPSETGENES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92979dc390984c7c638eadfba310ea85a821e56bb39fffb7833c872688eb1532/data'
        : 'community.wave.seqera.io/library/doublet_removal:d3212f3a1d58d039'}"

    input:
    tuple val(meta), val(names), path(h5ads)

    output:
    tuple val(meta), path("*.png"), emit: plot, optional: true
    path ("*_mqc.json"), emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

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
