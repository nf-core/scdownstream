process CUSTOM_COLLECTSIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas:2.2.2--4f1ae4025b6a600d':
        'community.wave.seqera.io/library/pandas:2.2.2--bd3db773995db54e' }"

    input:
    tuple val(meta), path(sizes)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path("*_mqc.json")            , emit: multiqc_files
    path "versions.yml"           , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'collectsizes.py'
}
