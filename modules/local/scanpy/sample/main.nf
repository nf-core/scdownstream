process SCANPY_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed34e104367da9d7fa40b39d78d86dce36948d95dbd2571abe923bc857bf192f/data' :
        'community.wave.seqera.io/library/python_pyyaml_scanpy:ec0559841ac88fb2' }"

    input:
    tuple val(meta), path(h5ad)
    val n
    val fraction

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"		   	   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_sampled"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('sample.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_sampled"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
