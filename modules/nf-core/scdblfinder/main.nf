// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SCDBLFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/709d3cc1733b4f3ddb44a51d7e52dafba6837b495886864bebc53690d95b5f6e/data'
        : 'community.wave.seqera.io/library/scdblfinder:03c14e340018af5c'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.predictions.parquet"), emit: predictions
    tuple val(meta), path("${prefix}.scores.parquet"), emit: scores
    path "versions.yml", emit: versions_scdblfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_scdblfinder"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('scdblfinder.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_scdblfinder"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.h5ad
    touch ${prefix}.predictions.parquet
    touch ${prefix}.scores.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        scdblfinder: \$(Rscript -e "suppressMessages(library(scDblFinder)); cat(as.character(packageVersion('scDblFinder')))")
        anndatar: \$(Rscript -e "suppressMessages(library(anndataR)); cat(as.character(packageVersion('anndataR')))")
        singlecellexperiment: \$(Rscript -e "suppressMessages(library(SingleCellExperiment)); cat(as.character(packageVersion('SingleCellExperiment')))")
        biocparallel: \$(Rscript -e "suppressMessages(library(BiocParallel)); cat(as.character(packageVersion('BiocParallel')))")
        rhdf5: \$(Rscript -e "cat(as.character(packageVersion('rhdf5')))")
        arrow: \$(Rscript -e "cat(as.character(packageVersion('arrow')))")
    END_VERSIONS
    """
}
