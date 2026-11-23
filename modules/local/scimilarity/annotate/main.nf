process SCIMILARITY_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4be3b0a53a34f772027d3132d99e71448cbfe0092cfcb034a45c9ab66cf30721/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/0648cc3712e639ff28425b32ce23ce9426fbb05c3bf0b6980a1cfe7e7ba2c23d/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/python_anndata_pyarrow_pyyaml_pruned:f4f0615a5198942f'
            : 'community.wave.seqera.io/library/python_anndata_pyarrow_pyyaml_pruned:bb1a128a5af6fc5c')}"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.parquet"), emit: obs
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    template('annotate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.parquet"
    touch "versions.yml"
    """
}
