process SCVITOOLS_SOLO {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/0703860e978cf8a1f420cd6f39a5d32a3e3bb255fbf547ae47f673ebc5b76fbe/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/4867a52f4a6f4a5216a3c21e2aead3f7eee65495b517c7b332e3bdf8cea13b30/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-gpu_cuda-version@sha256:069836e30ad5879e2ca3d2fc1c0904fdbff649add0cd19ac3c31afddc55aee4e'
            : 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-cpu:f0c65d8a2e3eacb6')}"

    input:
    tuple val(meta), path(h5ad)
    val batch_key
    val max_epochs

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.parquet"), emit: predictions
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('solo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=./tmp

    touch ${prefix}.h5ad
    touch ${prefix}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvi: \$(python3 -c 'import scvi; print(scvi.__version__)')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """
}
