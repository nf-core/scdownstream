process CUSTOM_COLLECTSIZES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0d347e0d943cb8d0624a8c73247db0c8e96e7c130fe0fb3125aeb7d1b140e79/data'
        : 'community.wave.seqera.io/library/anndata_python_pyyaml:5ac0fd4d280528a5'}"

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
