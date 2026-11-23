process ADATA_READCSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/378d6d9781b5231390b501a6c1c89576dd3ca8c19f6634a0325793c05ded8901/data'
        : 'community.wave.seqera.io/library/anndata_python_pyyaml:f7811e63d92e9821'}"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('readcsv.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
