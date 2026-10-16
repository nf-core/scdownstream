#!/usr/bin/env python3

import os
import re
import runpy
import sys
import tempfile
import unittest
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]


def dense(matrix):
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def fixture(name, matrix=None):
    if matrix is None:
        matrix = np.random.default_rng(42).poisson(3, size=(24, 30)).astype(np.float32)
    return ad.AnnData(
        sparse.csr_matrix(matrix),
        obs=pd.DataFrame(
            {"batch": name, "label": "Unknown"},
            index=[f"{name}_{i}" for i in range(len(matrix))],
        ),
        var=pd.DataFrame(index=[f"G{i}" for i in range(matrix.shape[1])]),
    )


def run_template(template, **variables):
    variables = {"task.process": "TEST", "task.cpus": "1", **variables}
    template_directory = os.environ.get("EXPRESSION_TEMPLATE_DIR")
    source_path = Path(template_directory) / Path(template).name if template_directory else ROOT / template
    source = source_path.read_text()
    source = re.sub(r"\$\{([^}]+)\}", lambda match: str(variables[match[1]]), source)
    script = Path("rendered.py")
    script.write_text(source)
    runpy.run_path(str(script), run_name="__main__")


def normalise(adata, prefix):
    source = f"{prefix}_input.h5ad"
    adata.write_h5ad(source)
    run_template("modules/local/scanpy/log1pnorm/templates/log1pnorm.py", h5ad=source, prefix=prefix)
    return ad.read_h5ad(f"{prefix}.h5ad")


def merge(reference, query):
    reference.write_h5ad("reference.h5ad")
    query.write_h5ad("query.h5ad")
    run_template(
        "modules/local/adata/merge/templates/merge.py",
        h5ads="query.h5ad",
        base="reference.h5ad",
        force_obs_cols="",
        prefix="merged",
    )
    return ad.read_h5ad("merged_outer.h5ad"), ad.read_h5ad("merged_integrate.h5ad")


class ExpressionLayersTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.previous_directory = Path.cwd()
        os.chdir(self.directory.name)
        self.addCleanup(self.directory.cleanup)
        self.addCleanup(os.chdir, self.previous_directory)

    def test_normalisation_preserves_counts_and_x(self):
        raw = np.array([[2.5, 0, 1], [5, 3, 2], [1, 2, 1]], dtype=np.float32)
        expected_log1p = np.log1p(raw * np.median(raw.sum(axis=1)) / raw.sum(axis=1)[:, None])
        for constructor in (np.asarray, sparse.csr_matrix, sparse.csc_matrix):
            for existing_counts in (False, True):
                with self.subTest(matrix=constructor.__name__, existing_counts=existing_counts):
                    adata = fixture("reference", raw)
                    adata.X = constructor(raw.copy())
                    if existing_counts:
                        adata.layers["counts"] = constructor(raw.copy())
                        adata.X = constructor(raw * 2)
                    original_x = dense(adata.X).copy()
                    result = normalise(adata, "normalised")
                    np.testing.assert_array_equal(dense(result.X), original_x)
                    np.testing.assert_array_equal(dense(result.layers["counts"]), raw)
                    np.testing.assert_allclose(dense(result.layers["log1p"]), expected_log1p, rtol=1e-6)

    def test_merge_layered_reference_and_layerless_query_preserves_counts(self):
        reference = fixture("reference")
        reference.layers["counts"] = reference.X.copy()
        reference.layers["log1p"] = reference.X.copy().log1p()
        query = fixture("query")
        query.var_names = [*query.var_names[:-1], "query_only"]
        outer, integrate = merge(reference, query)
        for result in (outer[query.obs_names, query.var_names], integrate[:, query.var_names]):
            np.testing.assert_array_equal(dense(result.layers["counts"]), dense(query.X))
        np.testing.assert_array_equal(
            dense(outer[reference.obs_names, reference.var_names].layers["counts"]), dense(reference.X)
        )

    def test_normalised_query_survives_feature_selection_and_reaches_scvi_as_counts(self):
        reference = normalise(fixture("reference"), "reference_normalised")
        query = fixture("query")
        query_raw = dense(query.X).copy()
        query = normalise(query, "query_normalised")
        outer, integrate = merge(reference, query)
        np.testing.assert_array_equal(dense(integrate.layers["counts"]), query_raw)
        np.testing.assert_allclose(dense(integrate.layers["log1p"]), dense(query.layers["log1p"]))
        np.testing.assert_array_equal(
            dense(outer[reference.obs_names].layers["log1p"]), dense(reference.layers["log1p"])
        )
        outer.write_h5ad("before_hvg.h5ad")
        run_template(
            "modules/local/scanpy/hvgs/templates/hvgs.py",
            h5ad="before_hvg.h5ad",
            prefix="hvg",
            n_hvgs=10,
            batch_key="",
            input_layer="log1p",
            log_normalize="false",
            excluded_genes="",
        )
        selected = ad.read_h5ad("hvg.h5ad")
        self.assertLess(selected.n_vars, outer.n_vars)
        expected = dense(outer[:, selected.var_names].X)
        np.testing.assert_array_equal(dense(selected.X), expected)
        self.assertTrue(np.all(dense(selected[query.obs_names].X).sum(axis=1) > 0))

        registered_matrices = []

        class CheckingSCVI:
            @staticmethod
            def setup_anndata(adata, **kwargs):
                registered = adata.layers[kwargs["layer"]] if kwargs.get("layer") else adata.X
                np.testing.assert_array_equal(dense(registered), expected)
                registered_matrices.append(registered)

            def __init__(self, adata, **kwargs):
                self.adata = adata

            def train(self, **kwargs):
                pass

            def get_latent_representation(self):
                return np.zeros((self.adata.n_obs, 2), dtype=np.float32)

            def save(self, path):
                pass

        scvi = ModuleType("scvi")
        scvi.settings = ModuleType("settings")
        scvi.__version__ = "test"
        model = ModuleType("scvi.model")
        model.SCVI = CheckingSCVI
        torch = ModuleType("torch")
        torch.set_float32_matmul_precision = lambda _: None
        torch.use_deterministic_algorithms = lambda _: None
        with patch.dict(sys.modules, {"scvi": scvi, "scvi.model": model, "torch": torch}):
            run_template(
                "modules/local/scvitools/scvi/templates/scvi.py",
                h5ad="hvg.h5ad",
                prefix="scvi",
                **{
                    "meta2.id ?: ''": "",
                    "categorical_covariates": "",
                    "continuous_covariates": "",
                    "batch_col": "batch",
                    "n_hidden": 8,
                    "n_layers": 1,
                    "n_latent": 2,
                    "dispersion": "gene",
                    "gene_likelihood": "nb",
                    "use_observed_lib_size": "true",
                    "task.accelerator ? 'true' : 'false'": "false",
                    "max_epochs": 1,
                    "max_epochs?:''": 1,
                },
            )
        self.assertEqual(len(registered_matrices), 1)


if __name__ == "__main__":
    unittest.main()
