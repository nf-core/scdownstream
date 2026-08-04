process SCANPY_PAGA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f36c00fcbaa1da1aaa6cf25cd8a562efb6728457d4a54b46fa06395aebd2083e/data' :
        'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:53120062012176bc' }"

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
