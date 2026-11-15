"""Fail-closed regression checks for the publication driver (no certificate rerun)."""
import copy
import hashlib
import json
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import verify_all as driver


class DriverChecks(unittest.TestCase):
    def setUp(self):
        self.profile = json.loads((driver.D4 / "data/census_floors.json").read_text())["profile"]

    def test_committed_group_domain(self):
        self.assertEqual(len(driver.census_groups(self.profile)), 124)

    def test_nonintegral_and_inexact_counts_fail(self):
        for value in ("1/2", 0.5, True, -1, 27):
            with self.subTest(value=value):
                altered = copy.deepcopy(self.profile)
                altered[0]["N"] = value
                with self.assertRaises(ValueError):
                    driver.census_groups(altered)

    def test_duplicate_missing_and_out_of_domain_fail(self):
        duplicate = copy.deepcopy(self.profile)
        duplicate[-1] = copy.deepcopy(duplicate[0])
        invalid = copy.deepcopy(self.profile)
        invalid[0].update(N=26, first_radius=1)
        for altered in (duplicate, self.profile[:-1], invalid):
            with self.assertRaises(ValueError):
                driver.census_groups(altered)

    def test_inexact_or_negative_floor_fails(self):
        for value in (0.25, True, "-1/100"):
            altered = copy.deepcopy(self.profile)
            altered[0]["new_floor"] = value
            with self.assertRaises(ValueError):
                driver.census_groups(altered)

    def test_optimized_interpreter_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "must-not-exist"
            run = subprocess.run([sys.executable, "-O", str(Path(driver.__file__)), "--out", str(out)],
                                 capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("disables assertions", run.stderr)
            self.assertFalse(out.exists())

    def test_export_manifests_are_current(self):
        self.assertGreater(driver.check_export_manifests(), 200)
        root = driver.D4 / "verify"
        for manifest_path in root.rglob("*MANIFEST.json"):
            manifest = json.loads(manifest_path.read_text())
            rows = manifest.get("files", [])
            if isinstance(rows, dict):
                rows = [dict(value, path=path) for path, value in rows.items()]
            for row in rows:
                path = manifest_path.parent / row["path"]
                with self.subTest(manifest=str(manifest_path), path=str(path)):
                    self.assertTrue(path.is_file())
                    data = path.read_bytes()
                    self.assertEqual(len(data), row["bytes"])
                    self.assertEqual(hashlib.sha256(data).hexdigest(), row["sha256"])

    def test_stage1_metadata_is_not_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "D4"
            shutil.copytree(driver.D4, root, ignore=shutil.ignore_patterns(".lake", "__pycache__", ".DS_Store"))
            stage = root / "verify/stage1_v7"
            spec = importlib.util.spec_from_file_location("metadata_fixture", stage / "run.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            (stage / "record/.DS_Store").write_bytes(b"OS metadata fixture")
            module.check_manifest()
            (stage / "record/unlisted.py").write_text("raise RuntimeError('must never run')\n")
            with self.assertRaises(module.Fail):
                module.check_manifest()


if __name__ == "__main__":
    unittest.main()
