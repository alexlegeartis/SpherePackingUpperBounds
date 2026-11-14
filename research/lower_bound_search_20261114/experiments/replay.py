"""Optional numerical experiments; these are NOT new minimum certificates."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import numpy as np

P = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("experiment", choices=["holes37", "new36", "section35", "p48_section38", "p48_bkz40_28",
                                        "KP34", "CTturyn1_scale_vertex", "from38_section37_0"])
args = parser.parse_args()
os.chdir(P)
subprocess.run(["g++", "-std=c++17", "-O3", "enum_shell.cpp", "-o", "enum_shell"], check=True)
name = "new36" if args.experiment == "holes37" else args.experiment
Q = np.load(name + ".npz")["Q"]
np.savetxt(name + "_R.txt", np.linalg.cholesky(Q).T)
bound = (160.000001 if name in ("section35", "new36") else
         8.00000001 if name in ("KP34", "CTturyn1_scale_vertex") else 6.00000001)
subprocess.run(["./enum_shell", str(len(Q)), str(bound), "8000000000", name + "_R.txt",
                name + "_vectors.txt", name + "_enum.log"], check=True)
if args.experiment == "holes37":
    subprocess.run([sys.executable, "holes37.py"], check=True)
else:
    subprocess.run([sys.executable, "eutaxy.py", name], check=True)
