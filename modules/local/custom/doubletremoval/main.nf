process CUSTOM_DOUBLETREMOVAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7f257d17ec38c213802f67cbe15cb9928b46d1aca4f492154b958f8b9195681c/data'
:         'community.wave.seqera.io/library/doublet_removal:be35baa3d64eb9b7' }"

    input:
    tuple val(meta), path(h5ad), path(predictions)
    val(threshold)
    val(removal)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*_mqc.json")             , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix    ?: "${meta.id}"
    template 'doubletremoval.py'

    stub:
    prefix = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
