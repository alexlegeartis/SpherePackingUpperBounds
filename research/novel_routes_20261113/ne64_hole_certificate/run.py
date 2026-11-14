"""Reproduce the complete Ne64 dual-hole interval certificate (four shards).

Run with Python from any directory. Requires NumPy, SymPy, and g++.
The recorded run used about 13 CPU minutes across four cores, after preparation.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import sys

P = Path(__file__).resolve().parent
SHARDS = 4
FLAGS = ["-O3", "-std=c++17", "-fno-fast-math", "-ffp-contract=off"]


def main():
    subprocess.run([sys.executable, str(P / "prepare.py")], check=True)
    subprocess.run(["g++", *FLAGS, str(P / "interval_enum.cpp"), "-o", str(P / "interval_enum")], check=True)
    names = ["setup.json", "ldl_bounds.txt", "gram.txt", "coset.txt", "interval_enum.cpp", "prepare.py"]
    manifest = {
        "shard_count": SHARDS,
        "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
        "compiler_flags": FLAGS,
        "input_and_source_sha256": {name: hashlib.sha256((P / name).read_bytes()).hexdigest() for name in names},
        "binary_sha256": hashlib.sha256((P / "interval_enum").read_bytes()).hexdigest(),
    }
    (P / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    def shard(i):
        args = [str(P / "interval_enum"), str(P / "ldl_bounds.txt"), str(P / "gram.txt"),
                str(P / "coset.txt"), str(i), str(SHARDS)]
        with (P / f"shard_{i}.json").open("w") as out, (P / f"shard_{i}.log").open("w") as log:
            subprocess.run(args, stdout=out, stderr=log, check=True)
        return i

    with ThreadPoolExecutor(max_workers=SHARDS) as pool:
        for i in pool.map(shard, range(SHARDS)):
            print(f"Shard {i} complete", flush=True)
    subprocess.run([sys.executable, str(P / "summarize.py")], check=True)
    subprocess.run([sys.executable, str(P.parent / "higher_parent" / "verify_ne65.py")], check=True)


if __name__ == "__main__":
    main()
