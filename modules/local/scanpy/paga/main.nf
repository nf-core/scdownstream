process SCANPY_PAGA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e2e626f2e43e75b0f218a91d16c13ac7efd15c6e0618dbc8ce06cd84218ca99/data'
        : 'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:1ea7adbfe9ea8cad'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path("*.pkl")                  , emit: uns, optional: true
    path("*.npy")                  , emit: obsp, optional: true
    path("*.png")                  , emit: plot, optional: true
    path("*_mqc.json")             , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    obs_key = meta.obs_key ?: "leiden"
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'paga.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "${prefix}_connectivities.npy"
    touch "${prefix}.png"
    touch "${prefix}_mqc.json"
    touch "versions.yml"
    """
}
