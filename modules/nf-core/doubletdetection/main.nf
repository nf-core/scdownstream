process DOUBLETDETECTION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/ddc1a22c6b4d0e7fdf1f855d1638f3952b6b88cebe2fde6dd59b97d300c4ae56/data'
        : 'community.wave.seqera.io/library/anndata_pip_pyarrow_python_pruned:af9720f5b3743d63'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.parquet"), emit: predictions
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('doubletdetection.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 -c 'import platform as pf; print(pf.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        doubletdetection: \$(python3 -c 'import doubletdetection as dt; print(dt.__version__)')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """
}
