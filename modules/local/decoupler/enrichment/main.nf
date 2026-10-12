process DECOUPLER_ENRICHMENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf16fb011f0dc4161a231eb90922f6fc850c70ab7679ccd01abb45357e4d29ff/data'
        : 'community.wave.seqera.io/library/decoupler_enrichment:5f69c22b99518995'}"

    input:
    tuple val(meta), path(de_results, stageAs: 'de/*')
    path networks, stageAs: 'networks/*'
    val method
    val tmin

    output:
    tuple val(meta), path("${prefix}_enrichment.parquet"), emit: results, optional: true
    path "${prefix}_decoupler_enrichment.pkl", emit: uns
    path "*_enrichment.png", emit: plots, optional: true
    path "*_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: (meta.de_method ? "${meta.id}_${meta.de_method}" : "${meta.id}")
    template('enrichment.py')

    stub:
    prefix = task.ext.prefix ?: (meta.de_method ? "${meta.id}_${meta.de_method}" : "${meta.id}")
    """
    touch ${prefix}_enrichment.parquet
    touch ${prefix}_decoupler_enrichment.pkl
    touch ${prefix}_network_enrichment.png
    touch ${prefix}_network_enrichment_mqc.json
    touch versions.yml
    """
}
