process DECOUPLER_ACTIVITY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5ca309507e02db5a619c02a9d8bb1d5b969ef6a815cceca6693784ee7fd061ec/data'
        : 'community.wave.seqera.io/library/decoupler_activity:8d8ed42929b0d74b'}"

    input:
    tuple val(meta), path(h5ad)
    path networks, stageAs: 'networks/*'
    val method
    val tmin

    output:
    tuple val(meta), path("${prefix}_*_scores.parquet"), emit: scores, optional: true
    path "${prefix}_decoupler_activity.pkl", emit: uns
    path "*_activity.png", emit: plots, optional: true
    path "*_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('activity.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_network_scores.parquet
    touch ${prefix}_decoupler_activity.pkl
    touch ${prefix}_network_activity.png
    touch ${prefix}_network_activity_mqc.json
    touch versions.yml
    """
}
