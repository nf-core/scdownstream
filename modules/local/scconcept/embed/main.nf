process SCCONCEPT_EMBED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a4f94db23f6a201f102b984db66975479051c43f493086cb2ad29bc487a335eb/data'
        : 'community.wave.seqera.io/library/build:b36ca665b4da8d71'}"

    input:
    tuple val(meta), path(h5ad, stageAs: 'input/counts.h5ad')
    tuple val(model_meta), path(config, stageAs: 'model/config.yaml'), path(checkpoint, stageAs: 'model/model.ckpt'), path(vocabulary, stageAs: 'model/gene_mappings/hsapiens.csv')

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("X_${prefix}.parquet"), emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('embed.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "X_${prefix}.parquet"
    touch "versions.yml"
    """
}
