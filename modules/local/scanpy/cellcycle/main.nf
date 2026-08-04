process SCANPY_CELLCYCLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed34e104367da9d7fa40b39d78d86dce36948d95dbd2571abe923bc857bf192f/data' :
        'community.wave.seqera.io/library/python_pyyaml_scanpy:ec0559841ac88fb2' }"

    input:
    tuple val(meta), path(h5ad)
    path s_genes
    path g2m_genes
    val symbol_col

    output:
    tuple val(meta), path("${prefix}.pkl"),  emit: obs
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "*_mqc.json"                        , emit: multiqc_files
    path "versions.yml",                     emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('cellcycle.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pkl
    touch ${prefix}.h5ad
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
