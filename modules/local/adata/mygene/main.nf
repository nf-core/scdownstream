process ADATA_MYGENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f121d9538b02d70d9d7770a37b50c34f7e6db07b57699b00877c8227927f6054/data' :
        'community.wave.seqera.io/library/mygene_anndata_pyyaml:2078fdff1608539a' }"

    input:
    tuple val(meta), path(h5ad)
    val(species)
    val(input_col)
    val(output_col)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

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
