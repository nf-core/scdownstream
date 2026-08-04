process ADATA_ENTROPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/354081f1aea522edf9eb2503417bfd3ac227e5f6cc419192042d47a1457bfb15/data' :
        'community.wave.seqera.io/library/python_pyyaml_anndata_scanpy:7df6225e45ce62d7' }"

    input:
    tuple val(meta), path(h5ad)
    val(group_col)
    val(entropy_col)
    val(plot_basis)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}.pkl"                   , emit: obs
    path "${prefix}.png"                   , emit: plots, optional: true
    path "${prefix}_mqc.json"              , emit: multiqc_files, optional: true
    path "versions.yml"                    , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    template 'entropy.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    if [ ${plot_basis ? 'true' : 'false'} ]; then
        touch ${prefix}.png
        touch ${prefix}_mqc.json
    fi

    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
