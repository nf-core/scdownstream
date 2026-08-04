process CUSTOM_DOUBLETREMOVAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1cda35bb187ad93a5ad29cd0a0340e1c4c5e3e84e8189dbbfca98e3bd4eef16b/data' :
        'community.wave.seqera.io/library/doublet_removal:81c2fac7b14cd21d' }"

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
