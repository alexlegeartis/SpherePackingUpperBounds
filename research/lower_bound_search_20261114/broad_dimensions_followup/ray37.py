"""Audit selected points on the 37D first-order ray by CVP upper witnesses."""

from pathlib import Path
import ctypes, json, math, subprocess
import numpy as np

P = Path(__file__).resolve().parent
work = P / "work"
work.mkdir(exist_ok=True)
subprocess.run(
    [
        "g++",
        "-O3",
        "-fPIC",
        "-shared",
        "-std=c++17",
        str(P / "cvp36.cpp"),
        "-o",
        str(work / "cvp36.so"),
    ],
    check=True,
)
z = np.load(P / "metric37_direction.npz")
Q = z["Q"]
L = z["L"]
D = z["D"]
E = L @ D @ L.T
G = Q[:36, :36]
R = np.ascontiguousarray(np.linalg.cholesky(G).T)
lib = ctypes.CDLL(str(work / "cvp36.so"))
pd = ctypes.POINTER(ctypes.c_double)
pi = ctypes.POINTER(ctypes.c_int64)
lib.init_cvp.argtypes = [ctypes.c_int, pd, pd]
lib.nearest_cvp.argtypes = [pd, ctypes.c_uint64, pd, pi]
lib.nearest_cvp.restype = ctypes.c_uint64
lib.init_cvp(36, G.copy().ctypes.data_as(pd), R.ctypes.data_as(pd))
out = []
limit = -1 / np.linalg.eigvalsh(D).min()
for fraction in [0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8]:
    t = limit * fraction
    Qn = Q + t * E
    Gn = np.ascontiguousarray(Qn[:36, :36])
    Rn = np.ascontiguousarray(np.linalg.cholesky(Gn).T)
    lib.init_cvp(36, Gn.ctypes.data_as(pd), Rn.ctypes.data_as(pd))
    s = np.linalg.solve(Gn, Qn[:36, 36])
    h2 = Qn[36, 36] - s @ Gn @ s
    minimum_upper = 8.0
    witnesses = []
    for j in [1, 2]:
        x = np.ascontiguousarray(j * s)
        v = ctypes.c_double()
        w = np.zeros(36, dtype=np.int64)
        nodes = lib.nearest_cvp(
            x.ctypes.data_as(pd), 500000, ctypes.byref(v), w.ctypes.data_as(pi)
        )
        minimum_upper = min(minimum_upper, v.value + j * j * h2)
        witnesses.append(
            {
                "layer": j,
                "distance_upper": v.value,
                "nearest_vector": w.tolist(),
                "nodes": int(nodes),
            }
        )
    det = math.exp(np.linalg.slogdet(Qn)[1])
    conditional = 2**18.5 / math.sqrt(det)
    upper = (minimum_upper / 4) ** 18.5 / math.sqrt(det)
    row = {
        "parameter": t,
        "conditional_center_density": conditional,
        "minimum_squared_upper": minimum_upper,
        "density_upper_from_witnesses": upper,
        "witnesses": witnesses,
    }
    out.append(row)
    print({k: v for k, v in row.items() if k != "witnesses"}, flush=True)
(P / "ray37_results.json").write_text(json.dumps(out, indent=2) + "\n")
