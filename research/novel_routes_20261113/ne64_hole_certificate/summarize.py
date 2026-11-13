"""Check that the exhaustive interval shards form a complete empty traversal.

The compact result is a computation record, not a substitute for rerunning
interval_enum.cpp when independently reproducing the exhaustive certificate.
"""
from pathlib import Path
import hashlib
import json

P = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    setup = json.loads((P / "setup.json").read_text())
    assert setup["inclusive_integer_norm_bound"] == 60
    assert setup["congruence_modulus"] == 3
    assert setup["residue_norm_mod_18"] == 6
    for name, digest in setup["sha256"].items():
        assert sha(P / name) == digest
    manifest = json.loads((P / "run_manifest.json").read_text())
    for name, digest in manifest["input_and_source_sha256"].items():
        assert sha(P / name) == digest
    n = manifest["shard_count"]
    shards = [json.loads((P / f"shard_{i}.json").read_text()) for i in range(n)]
    for i, row in enumerate(shards):
        assert row["shard_id"] == i and row["shard_count"] == n
        assert row["status"] == "COMPLETE" and row["found"] is False
        assert row["split_depth"] == 47
        assert row["norm"] == -1 and row["witness"] == []
        assert row["nodes"] > 0
        assert row["assigned_prefixes"] == (row["prefixes"] + n - 1 - i) // n
    counts = {row["prefixes"] for row in shards}
    assert len(counts) == 1
    prefixes = counts.pop()
    assert sum(row["assigned_prefixes"] for row in shards) == prefixes
    out = {
        "status": "complete_outward_interval_enumeration",
        "certified_empty": True,
        "dimension": 64,
        "class_id": 4,
        "inclusive_integer_norm_bound": 60,
        "congruence_modulus": 3,
        "residue_norm_mod_18": 6,
        "certified_coset_minimum_lower_bound": "26/3",
        "method": "Exact rational LDL enclosed in dyadic intervals; exhaustive IEEE-754 outward-rounded recursion; exact integer leaves",
        "arithmetic_assumptions": [
            "IEEE-754 binary64 round-to-nearest, with nextafter outward bounds",
            "Compilation without unsafe fast-math or floating-point contraction",
            "Signed 128-bit integer operations stay within the checked input and coordinate caps",
        ],
        "shard_count": n,
        "all_shards_complete": True,
        "total_nodes_including_repeated_upper_tree": sum(row["nodes"] for row in shards),
        "total_exact_leaves": sum(row["leaves"] for row in shards),
        "shared_prefix_count": prefixes,
        "sum_assigned_prefixes": sum(row["assigned_prefixes"] for row in shards),
        "shards": shards,
        "input_sha256": setup["sha256"],
        "proof_source_sha256": sha(P / "interval_enum.cpp"),
        "prepare_source_sha256": sha(P / "prepare.py"),
        "run_manifest_sha256": sha(P / "run_manifest.json"),
        "scope": "The coset computation is exhaustive. The parent Ne64 minimum 12 is inherited from the published theorem. Global sphere-packing priority is not asserted.",
    }
    (P / "summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("shards", "input_sha256")}, indent=2))


if __name__ == "__main__":
    main()
