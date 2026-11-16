process SCVITOOLS_SOLO {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/0072fc1c3e82f329a0f92c8544e6f9b2129ff08f7cbfa32b50f8c3df1b54f3f3/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba88d310e8f5d69d321b94bb73012e8f5571ec4aee16e1047e7dd378b606fb2a/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-gpu_cuda-version@sha256:7bd8b3dd714d9e419a67bbbaab1d92fa5823fd370b307c7e67d38a02ed36cbb1'
            : 'community.wave.seqera.io/library/pyarrow_scvi-tools:002a45f506b3d900')}"

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
