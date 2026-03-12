process PHATE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://miguelrosell/nf-core-manifold:dev' :
        'docker.io/miguelrosell/nf-core-manifold:dev' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ./tmp
    export MPLCONFIGDIR="./tmp"
    export NUMBA_CACHE_DIR="./tmp"

    run_phate.py \\
        --input ${h5ad} \\
        --output ${prefix}_phate.h5ad \\
        --n_jobs ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phate: \$(python -c "import importlib.metadata; print(importlib.metadata.version('phate'))")
        scanpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('scanpy'))")
    END_VERSIONS
    """
}
