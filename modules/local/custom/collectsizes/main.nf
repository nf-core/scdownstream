process CUSTOM_COLLECTSIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3ff2a833e55ecf43bf445ce4a63faa9388ee7fceb8d51e603479572694b6b9b4/data' :
        'community.wave.seqera.io/library/custom_collectsizes:5b6a03eb50fdf124' }"

    input:
    tuple val(meta), path(sizes)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path("*_mqc.json")            , emit: multiqc_files
    path "versions.yml"           , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'collectsizes.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sizes.tsv
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
