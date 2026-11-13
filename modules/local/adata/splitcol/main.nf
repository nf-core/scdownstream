process ADATA_SPLITCOL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de2455bdcbc53723a8d942ffb0dec7559690ef4bea6b380fea2edba7c5d637b9/data'
:         'community.wave.seqera.io/library/anndata_python_pyyaml:e77bb4f393d60a08' }"

    input:
    tuple val(meta), path(h5ad)
    val column

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('split_column.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch A.h5ad
    touch B.h5ad
    touch versions.yml
    """
}
