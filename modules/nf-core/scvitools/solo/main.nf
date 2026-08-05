process SCVITOOLS_SOLO {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/21c048694d2a9bc03c17a310c79c5691d2fc83df62e5e2654791a2105da8b09f/data'
:         'community.wave.seqera.io/library/scvi-tools_pyarrow:1734dedd3c3d134b' }"

    input:
    tuple val(meta), path(h5ad)
    val(batch_key)
    val(max_epochs)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.parquet") , emit: predictions
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'solo.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=./tmp

    touch ${prefix}.h5ad
    touch ${prefix}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvi: \$(python3 -c 'import scvi; print(scvi.__version__)')
    END_VERSIONS
    """
}
