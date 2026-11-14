process CUSTOM_DOUBLETREMOVAL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56e5b2e0a7deae4516efa83ada0d514260d70199270fdd7d66e0ed7841bc3f27/data'
        : 'community.wave.seqera.io/library/doublet_removal:3e0acc996063ca2e'}"

    input:
    tuple val(meta), path(h5ad), val(methods), path(predictions, stageAs: 'predictions/*', arity: '1..*')
    val threshold
    val removal

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.obs.pkl"), emit: obs
    path ("*_mqc.json"), emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('doubletremoval.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.obs.pkl
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
