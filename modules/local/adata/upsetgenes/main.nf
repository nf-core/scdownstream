process ADATA_UPSETGENES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45339bf761a2cf0cdb058492bc37f3df8b05b363731d491d1d3a14e9ba0b8f55/data'
        : 'community.wave.seqera.io/library/harmonypy_anndata_leidenalg_numpy_pruned:43066d5f86f18261'}"

    input:
    tuple val(meta), val(names), path(h5ads)

    output:
    tuple val(meta), path("*.png"), emit: plot, optional: true
    path ("*_mqc.json")           , emit: multiqc_files, optional: true
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('upsetplot.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
