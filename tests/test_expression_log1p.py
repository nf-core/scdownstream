"""Run in the rankgenesgroups environment with python -m unittest discover -s tests -p test_expression_log1p.py."""

import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse


class BoundaryReachedError(Exception):
    pass


class TestExpressionLog1p(unittest.TestCase):
    def test_expression_at_analysis_or_export_boundary(self):
        root = Path(__file__).resolve().parents[1] / "modules/local"
        counts = np.tile(np.array([[0, 1, 8], [3, 0, 15]], dtype=np.float32), (4, 1))
        for module, template, target in (
            ("scanpy/rankgenesgroups", "rank_genes_groups.py", "scanpy.tl.rank_genes_groups"),
            ("adata/prepcellxgene", "prepcellxgene.py", "anndata.AnnData.write_h5ad"),
        ):
            for sparse in (False, True):
                for apply_log1p, logged in ((True, False), (False, True), (False, False)):
                    with self.subTest(module=module, sparse=sparse, apply_log1p=apply_log1p, logged=logged):
                        values = np.log1p(counts) if logged else counts.copy()
                        expected = np.log1p(values) if apply_log1p else values.copy()
                        adata = ad.AnnData(
                            X=csr_matrix(values) if sparse else values,
                            obs=pd.DataFrame({"group": ["a", "b"] * 4}, index=[f"cell_{i}" for i in range(8)]),
                        )
                        if logged:
                            adata.uns["log1p"] = {"base": None}
                        captured = []

                        def capture(actual, *args, **kwargs):
                            matrix = actual.X.toarray() if issparse(actual.X) else actual.X
                            captured.append(matrix.copy())
                            raise BoundaryReachedError

                        with tempfile.TemporaryDirectory() as directory:
                            h5ad = Path(directory) / "input.h5ad"
                            adata.write_h5ad(h5ad, compression="gzip")
                            substitutions = {
                                "h5ad": str(h5ad),
                                "prefix": str(Path(directory) / "result"),
                                "obs_key": "group",
                                "method": "wilcoxon",
                                "rank_key": "rank_genes_groups",
                                "filter_col ?: ''": "",
                                "filter_val ?: ''": "",
                                "meta.id": "test",
                                "meta.integration": "test",
                                "task.cpus": "1",
                                "task.process": "EXPRESSION_TEST",
                                "apply_log1p": str(apply_log1p).lower(),
                            }
                            source = (root / module / "templates" / template).read_text()
                            source = re.sub(r"\$\{([^}]+)\}", lambda match: substitutions[match[1]], source)
                            with (
                                patch(target, autospec=True, side_effect=capture),
                                self.assertRaises(BoundaryReachedError),
                            ):
                                exec(compile(source, template, "exec"), {"__name__": "__main__"})
                        self.assertEqual(len(captured), 1)
                        np.testing.assert_allclose(captured[0], expected, rtol=1e-6, atol=0)


if __name__ == "__main__":
    unittest.main()
