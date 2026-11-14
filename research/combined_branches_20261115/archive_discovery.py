"""Archive complete continuation history after SHA256 read-back verification."""

from pathlib import Path
import argparse, hashlib, json, tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--output", default="discovery_history.tar.xz")
ap.add_argument("--manifest", default="DISCOVERY_MANIFEST.json")
ap.add_argument(
    "--dimensions", type=int, nargs="+", default=[26, 27, 28, 29, 30, 31, 33, 34]
)
args = ap.parse_args()
archive = HERE / args.output
assert not archive.exists(), "Choose a new archive filename"
retained = []
members = []
keep = {
    "candidate.json",
    "angular.json",
    "angular_verified.json",
    "radial_certificate.json",
    "radial_export.json",
    "cuts.json",
}
for n in args.dimensions:
    for directory in (ROOT / "research/ambitious_profiles_20261114").glob(
        f"integrated_n{n}_*"
    ):
        for file in sorted(directory.rglob("*")):
            if not file.is_file():
                continue
            assert not file.is_symlink()
            row = {
                "path": str(file.relative_to(ROOT)),
                "sha256": hashlib.sha256(file.read_bytes()).hexdigest(),
                "bytes": file.stat().st_size,
            }
            if file.parent == directory and file.name in keep:
                retained.append(row)
            else:
                members.append(row)
with tarfile.open(archive, "w:xz", preset=6) as out:
    for row in members:
        out.add(ROOT / row["path"], arcname=row["path"], recursive=False)
with tarfile.open(archive, "r:xz") as inp:
    actual = inp.getmembers()
    assert len(actual) == len(members)
    for row, entry in zip(members, actual):
        assert entry.isfile() and entry.name == row["path"]
        assert (
            hashlib.sha256(inp.extractfile(entry).read()).hexdigest() == row["sha256"]
        )
report = {
    "archive": str(archive.relative_to(ROOT)),
    "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    "archive_members": members,
    "retained_latest": retained,
    "all_members_read_back_verified": True,
    "extraction_directory": "repository root",
}
(HERE / args.manifest).write_text(json.dumps(report, indent=2) + "\n")
for row in members:
    (ROOT / row["path"]).unlink()
print(
    "PASS:",
    len(members),
    "checkpoint files archived and verified;",
    len(retained),
    "latest inputs retained",
)
