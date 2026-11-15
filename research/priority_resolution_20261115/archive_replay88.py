"""Bind a complete, freshly executed upstream 88D replay to our catalogue input.

This archives evidence of an upstream computation, not an independent
reimplementation or a standalone minimum verifier. See REPLAY88.md.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np

PIN = "628ce2a7ac23acf8d4c295677e4310c9ec5b21d9"
P = Path(__file__).resolve().parent


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", type=Path, required=True)
    ap.add_argument("--run", type=Path, required=True)
    args = ap.parse_args()
    u, r = args.upstream.resolve(), args.run.resolve()
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=u, text=True).strip() == PIN
    assert subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=u) == b""
    read = lambda name: json.loads((r/name).read_text())
    result = read("verification.json")
    assert result["status"] == "fully_reproduced"
    assert len(result["stages"]) == 8
    assert all(s["returncode"] == 0 for s in result["stages"])
    assert read("L/REPRODUCTION.json")["status"] == "fully_reproduced"
    assert read("L/FINAL_VERIFICATION.json")["status"] == "verified"
    assert read("L/FINAL_VERIFICATION.json")["trace_minimum"] == 8
    j = read("J0/output/rank4_audit/FULL_RUN_SUMMARY.json")
    assert j["verified"] and not j["rank4_h3_exists"]
    assert j["consequence_using_existing_rank1_to3_bounds_and_norm4_witness"]["trace_minimum"] == 8
    assert read("groups/verification.json")["verified"]
    for name in ("factor_claims", "L_perfection", "L0_perfection"):
        assert read(name+"/verification.json")["status"] == "verified"
    assert read("L0_input/lattice.json") == json.loads((u/"data/L0/lattice.json").read_text())
    assert read("L/workspace/output/perfect_reaudit/input/lattice.json") == json.loads((u/"data/L/lattice.json").read_text())
    g = read("L0_input/lattice.json")["G"]
    for n in (86, 87):
        with np.load(P.parent/f"cross_dimension_lower_20261115/certificates/claimed88_n{n}.npz", allow_pickle=False) as d:
            assert d["G"].tolist() == g
    branches = sorted((r/"L/workspace/output/d4_projective_J/full_run/anchors").glob("*/independent_certificate/verification.json"))
    assert len(branches) == 78
    assert all(json.loads(p.read_text())["status"] == "verified" for p in branches)
    target = P/"replay88"
    target.mkdir(exist_ok=True)
    files = ["verification.json", "stages.json", "L/REPRODUCTION.json", "L/FINAL_VERIFICATION.json",
             "J0/output/rank4_audit/FULL_RUN_SUMMARY.json", "groups/verification.json",
             "factor_claims/verification.json", "L_perfection/verification.json", "L0_perfection/verification.json"]
    for name in files:
        dest = target/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(r/name, dest)
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=u).decode().split("\0")
    manifest = {"source_repository": "https://github.com/RyotaroSakamoto0113/rank-88-extremal-lattice-computations",
                "commit": PIN, "status": "full_upstream_computational_proof_freshly_reproduced",
                "new_implementation": False, "formal_proof": False,
                "wall_seconds": result["wall_seconds"],
                "catalogue_parent": "ideal115L0", "catalogue_gram_equals_rebuilt_L0": True,
                "catalogue_parent_minimum_squared": 8,
                "upstream_tracked_file_sha256": {name: digest(u/name) for name in tracked if name},
                "fresh_report_sha256": {name: digest(r/name) for name in files},
                "fresh_lattice_sha256": {name: digest(r/name) for name in
                    ("L0_input/lattice.json", "L/workspace/output/perfect_reaudit/input/lattice.json")},
                "independent_branch_report_sha256": {str(p.relative_to(r)): digest(p) for p in branches}}
    (target/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("Full upstream replay bound to both archived antipode certificates; 78 branches verified.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O.")
    main()
