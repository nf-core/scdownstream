process SCANPY_RANKGENESGROUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3bfaa458a82be5af0d36e828f85c953789fa66f088bfc57b48e6f2fd6e77a79a/data' :
        'community.wave.seqera.io/library/adjusttext_python_pyyaml_scanpy:bb9fae9f02f816e2' }"

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    tuple val(filter_col), val(filter_val)
    val(method)
    val(rank_key)
    path interesting_genes

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path "*.pkl"                   , emit: uns, optional: true
    path "*.png"                   , emit: plots, optional: true
    path "*.csv"                   , emit: markers, optional: true
    path "*_mqc.json"              , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('rank_genes_groups.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "${prefix}.png"
    touch "${prefix}_dotplot.png"
    touch "${prefix}_volcano.png"
    touch "${prefix}_mqc.json"
    touch "${prefix}_dotplot_mqc.json"
    touch "${prefix}_volcano_mqc.json"
    touch "versions.yml"
    """
}
