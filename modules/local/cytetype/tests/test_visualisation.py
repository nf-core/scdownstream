"""Exercise the template's CyteType initialisation and submission boundary offline."""

import ast
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock


class VisualisationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        template = Path(__file__).resolve().parents[1] / "templates/cytetype.py"
        tree = ast.parse(template.read_text())
        start = next(
            i
            for i, node in enumerate(tree.body)
            if isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "annotator" for target in node.targets)
        )
        end = next(
            i
            for i, node in enumerate(tree.body[start:], start)
            if isinstance(node, ast.Assign)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Attribute)
            and isinstance(node.value.func.value, ast.Name)
            and node.value.func.value.id == "annotator"
            and node.value.func.attr == "run"
        )
        cls.submission = compile(ast.Module(body=tree.body[start : end + 1], type_ignores=[]), str(template), "exec")

    def submit(self, visualisation):
        def run(**kwargs):
            self.assertEqual(kwargs["study_context"], "Human PBMCs")
            self.assertIsInstance(annotator.visualization_data["coordinates"], list)
            self.assertIsInstance(annotator.visualization_data["clusters"], list)
            self.assertEqual(
                len(annotator.visualization_data["coordinates"]), len(annotator.visualization_data["clusters"])
            )
            return adata

        adata = object()
        annotator = SimpleNamespace(visualization_data=visualisation, run=Mock(side_effect=run))
        factory = Mock(return_value=annotator)
        namespace = {
            "CyteType": factory,
            "adata_work": adata,
            "group_key": "leiden",
            "rank_key": "rank_genes_groups",
            "auth_token_arg": None,
            "study_context": "Human PBMCs",
        }
        exec(self.submission, namespace)
        factory.assert_called_once_with(adata, group_key="leiden", rank_key="rank_genes_groups", auth_token=None)
        annotator.run.assert_called_once()
        self.assertIs(namespace["adata_work"], adata)
        return annotator.visualization_data

    def test_missing_coordinates_use_empty_arrays(self):
        result = self.submit({"coordinates": None, "clusters": []})
        self.assertEqual(result, {"coordinates": [], "clusters": []})

    def test_real_coordinates_and_cluster_labels_are_preserved(self):
        visualisation = {"coordinates": [[1.0, 2.0], [3.0, 4.0]], "clusters": ["0", "1"]}
        result = self.submit(visualisation)
        self.assertIs(result, visualisation)
        self.assertEqual(result, {"coordinates": [[1.0, 2.0], [3.0, 4.0]], "clusters": ["0", "1"]})

    def test_empty_visualisation_is_preserved(self):
        visualisation = {"coordinates": [], "clusters": []}
        self.assertIs(self.submit(visualisation), visualisation)


if __name__ == "__main__":
    unittest.main()
