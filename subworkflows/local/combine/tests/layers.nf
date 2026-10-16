include { COMBINE           } from '../main'
include { FEATURE_SELECTION } from '../../feature_selection'

process CREATE_EXPRESSION_FIXTURES {
    conda "${moduleDir}/../../../../modules/local/scanpy/log1pnorm/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    tuple path(regression), path(templates)

    output:
    path 'reference.h5ad', emit: reference
    path 'query*.h5ad', emit: queries

    script:
    """
    EXPRESSION_TEMPLATE_DIR="\$PWD" python ${regression} -v
    python - <<'PY'
    import sys
    sys.path.insert(0, '.')
    from expression_layers import dense, fixture
    import numpy as np

    reference = fixture('reference')
    reference.layers['counts'] = reference.X.copy()
    raw = dense(reference.X)
    reference.layers['log1p'] = np.log1p(raw * 10000 / raw.sum(axis=1)[:, None])
    reference.write_h5ad('reference.h5ad')
    fixture('query1').write_h5ad('query1.h5ad')
    query2 = fixture('query2')
    query2.X *= 3
    query2.var_names = [*query2.var_names[:-1], 'query_only']
    query2.write_h5ad('query2.h5ad')
    PY
    """
}

process CHECK_EXPRESSION_LAYERS {
    conda "${moduleDir}/../../../../modules/local/scanpy/log1pnorm/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    path regression
    path reference
    path queries
    tuple val(meta), path(outer)
    tuple val(meta2), path(selected)

    output:
    path 'passed.txt'

    script:
    """
    python - <<'PY'
    import sys
    sys.path.insert(0, '.')
    from expression_layers import dense
    import anndata as ad
    import numpy as np
    from pathlib import Path

    outer = ad.read_h5ad('${outer}')
    selected = ad.read_h5ad('${selected}')
    reference = ad.read_h5ad('${reference}')
    np.testing.assert_array_equal(dense(outer.layers['counts']), dense(outer.X))
    np.testing.assert_array_equal(dense(selected.layers['counts']), dense(selected.X))
    np.testing.assert_array_equal(
        dense(outer[reference.obs_names, reference.var_names].layers['log1p']),
        dense(reference.layers['log1p']),
    )
    for path in '${queries}'.split():
        query = ad.read_h5ad(path)
        raw = dense(query.X)
        expected_log1p = np.log1p(raw * np.median(raw.sum(axis=1)) / raw.sum(axis=1)[:, None])
        merged_query = outer[query.obs_names, query.var_names]
        np.testing.assert_array_equal(dense(merged_query.X), raw)
        np.testing.assert_array_equal(dense(merged_query.layers['counts']), raw)
        np.testing.assert_allclose(dense(merged_query.layers['log1p']), expected_log1p, rtol=1e-6)
        assert np.all(dense(selected[query.obs_names].X).sum(axis=1) > 0)
    assert selected.n_obs == outer.n_obs
    assert selected.n_vars < outer.n_vars
    Path('passed.txt').write_text('Expression layers preserved\\n')
    PY
    """
}

workflow TEST_EXTENSION_LAYERS {
    take:
    ch_sources

    main:
    CREATE_EXPRESSION_FIXTURES(ch_sources)
    ch_queries = CREATE_EXPRESSION_FIXTURES.out.queries
        .flatten()
        .map { h5ad -> [[id: h5ad.simpleName], h5ad] }

    COMBINE(
        ch_queries,
        CREATE_EXPRESSION_FIXTURES.out.reference.map { h5ad -> [[id: 'reference'], h5ad] },
        true,
        'hvgs',
        25,
        '',
        null,
        'log1p',
        null,
        null,
        '',
        '',
        8,
        1,
        2,
        'gene',
        'nb',
        1,
        true,
        null,
        null,
        null,
        '',
        false,
        0,
        'stratified_label_batch',
        42,
        'fast',
        '',
    )

    FEATURE_SELECTION(COMBINE.out.h5ad_inner, 'hvgs', 25, [], 'log1p')

    CHECK_EXPRESSION_LAYERS(
        ch_sources.map { regression, _templates -> regression },
        CREATE_EXPRESSION_FIXTURES.out.reference,
        CREATE_EXPRESSION_FIXTURES.out.queries,
        COMBINE.out.h5ad,
        FEATURE_SELECTION.out.h5ad,
    )

    emit:
    CHECK_EXPRESSION_LAYERS.out
}
