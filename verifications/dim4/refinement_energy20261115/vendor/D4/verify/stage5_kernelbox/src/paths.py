"""Input and output locations for stage 5. Path plumbing only.

Every path is relative to this file, so the package runs from any folder.
Environment overrides (set by run.py):
  D4_STAGE2_DIR  folder holding ACCEPTANCE.json (default: verify/stage2_origin/reference)
  D4_CENSUS_FILE the census floors file (default: data/census_floors.json)
  D4_STAGE5_OUT  folder for auxiliary outputs (default: stage5_kernelbox/out)
"""
from pathlib import Path
import os

SRC = Path(__file__).resolve().parent
STAGE = SRC.parent
PKG = STAGE.parents[1]
DATA = PKG / "data"
V7 = DATA / "v7"
KBOX = DATA / "kernel-box"
STAGE2 = Path(os.environ["D4_STAGE2_DIR"]) if os.environ.get("D4_STAGE2_DIR") else PKG / "verify" / "stage2_origin" / "reference"
ACCEPTANCE = STAGE2 / "ACCEPTANCE.json"
CENSUS = Path(os.environ["D4_CENSUS_FILE"]) if os.environ.get("D4_CENSUS_FILE") else DATA / "census_floors.json"
OUT = Path(os.environ["D4_STAGE5_OUT"]) if os.environ.get("D4_STAGE5_OUT") else STAGE / "out"
