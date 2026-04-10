process ADATA_MYGENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95d94db7b6389558ec95d6ad737fd962d8ade08a2323f298c92afed4c240c177/data':
        'community.wave.seqera.io/library/mygene_anndata_pyyaml:d9454f09fb1f98d5' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    input_col = task.ext.input_col ?: "index"
    output_col = task.ext.output_col ?: "symbols"

    if ("${prefix}.h5ad" == "${h5ad}")
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    template 'mygene.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
