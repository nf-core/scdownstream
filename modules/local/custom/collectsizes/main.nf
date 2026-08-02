process CUSTOM_COLLECTSIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68c5022a4fc1d7a52acf4d7ce0201540ff4f72b4eea2c5598ab31326b98a1206/data' :
        'community.wave.seqera.io/library/custom_collectsizes:ddf7b396f2bfaa5d' }"

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
