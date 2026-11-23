process CUSTOM_COLLECTSIZES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/378d6d9781b5231390b501a6c1c89576dd3ca8c19f6634a0325793c05ded8901/data'
        : 'community.wave.seqera.io/library/anndata_python_pyyaml:f7811e63d92e9821'}"

    input:
    tuple val(meta), path(sizes)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path ("*_mqc.json"), emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('collectsizes.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sizes.tsv
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
