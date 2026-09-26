process ADATA_MYGENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/ccc5af513e241deec87c61698bd3d49f458a78938474ff6e2b9a006dfc51690a/data'
:         'community.wave.seqera.io/library/mygene_anndata_python_pyyaml:0b27d9b10dbf88bf' }"

    input:
    tuple val(meta), path(h5ad)
    val(species)

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
