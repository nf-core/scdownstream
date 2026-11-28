process DECOUPLER_ACTIVITY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2821b036063d6024ae937e9f6cbc1b777333ac3bdcca793f68e1fbe868ae976/data'
        : 'community.wave.seqera.io/library/decoupler_activity:374661b94d4383e8'}"

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
