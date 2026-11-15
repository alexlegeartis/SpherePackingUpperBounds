import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "verify_projected_eight", HERE / "src/verify_projected_eight.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


class ProjectedEightChecks(unittest.TestCase):
    def setUp(self):
        self.path = HERE.parents[1] / "data/projected_eight.json"
        self.data = json.loads(self.path.read_text())

    def test_exact_certificate_passes(self):
        result = CHECKER.verify(self.path)
        self.assertEqual(result["status"], "EXACT_PROJECTED_EIGHT_PASS")
        self.assertEqual(result["count_bound"], 8)
        self.assertEqual(result["objective"], "89804130521/9999900000")
        self.assertEqual(result["bernstein_leaves"], 12)
        self.assertEqual(result["bernstein_depth"], 7)

    def test_bad_cutoff_and_negative_coefficient_fail(self):
        for altered in (copy.deepcopy(self.data), copy.deepcopy(self.data)):
            if altered["cutoff"] == "223/735":
                altered["cutoff"] = "1/3"
            else:
                altered["coefficients"][1] = "-1/10000000000"
            with tempfile.NamedTemporaryFile("w", suffix=".json") as handle:
                json.dump(altered, handle)
                handle.flush()
                with self.assertRaises(ValueError):
                    CHECKER.verify(Path(handle.name))

    def test_float_coefficients_are_rejected(self):
        altered = copy.deepcopy(self.data)
        altered["coefficients"][1] = 0.1
        with tempfile.NamedTemporaryFile("w", suffix=".json") as handle:
            json.dump(altered, handle)
            handle.flush()
            with self.assertRaises(ValueError):
                CHECKER.verify(Path(handle.name))


if __name__ == "__main__":
    unittest.main()
