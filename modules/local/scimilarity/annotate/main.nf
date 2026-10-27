process SCIMILARITY_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c90e92bde5ea63596f8af5f32ca881d0dc1053bd7c8d1e484a70b1d3d6fef56d/data'
        : 'community.wave.seqera.io/library/python_anndata_pyyaml_zarr_pruned:726de88bc84fbfa1'}"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.pkl")                 , emit: obs
    path "versions.yml"                    , emit: versions, topic: versions

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
    touch "${prefix}.pkl"
    touch "versions.yml"
    """
}
