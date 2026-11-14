"""Floating-discovery catalog seed screen; source minimum is a literature input."""

from pathlib import Path
import numpy as np, re, html, json, subprocess, time
from flint import fmpz_mat
from fpylll import IntegerMatrix, GSO, LLL

P = Path(__file__).resolve().parent
O = P / "work" / "catalog"
O.mkdir(parents=True, exist_ok=True)
binary = O / "enum_shell"
if not binary.exists():
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++17",
            str(P.parent / "experiments" / "enum_shell.cpp"),
            "-o",
            str(binary),
        ],
        check=True,
    )


def field(s, key):
    match = re.search(
        r'<a\s+NAME="' + key + r'[^\"]*"[^>]*>.*?</a>\s*<br\s*/?>', s, re.S | re.I
    )
    if not match:
        return ""
    block = re.split(r"<a\s+NAME=", s[match.end() :], flags=re.I)[0]
    return html.unescape(re.sub("<[^>]+>", " ", block)).strip()


rows = []
start = time.time()
for path in sorted((P / "catalog_sources").glob("*.html")):
    name = path.stem
    s = path.read_text()
    ns = re.findall(r"\b\d+\b", field(s, "DIMENSION"))
    ms = re.findall(r"\b\d+\b", field(s, "MINIMAL_NORM"))
    if not ns or not ms:
        continue
    n = int(ns[0])
    minimum = int(ms[0])
    raw = field(s, "GRAM")
    values = list(map(int, re.findall(r"-?\d+", raw)))
    G = None
    for offset in (2, 1, 0):
        if len(values) >= offset + n * n:
            A = np.array(values[offset : offset + n * n], dtype=np.int64).reshape(n, n)
            if np.array_equal(A, A.T) and np.linalg.eigvalsh(A.astype(float)).min() > 0:
                G = A
                break
        if len(values) == offset + n * (n + 1) // 2:
            A = np.zeros((n, n), dtype=np.int64)
            A[np.tril_indices(n)] = values[offset:]
            A = A + A.T - np.diag(A.diagonal())
            if np.linalg.eigvalsh(A.astype(float)).min() > 0:
                G = A
                break
    if G is None:
        continue
    det = int(fmpz_mat(G.tolist()).det())
    b = IntegerMatrix.from_matrix(G.tolist())
    u = IntegerMatrix.identity(n)
    gso = GSO.Mat(b, U=u, gram=True, flags=GSO.INT_GRAM)
    LLL.Reduction(gso)()
    U = np.array([[int(u[i, j]) for j in range(n)] for i in range(n)], dtype=np.int64)
    Q = U @ G @ U.T
    assert np.array_equal(
        Q, np.array([[int(b[i, j]) for j in range(n)] for i in range(n)])
    )
    prefix = O / name
    np.savez_compressed(str(prefix) + ".npz", Q=Q, U=U, original=G, minimum=minimum)
    np.savetxt(str(prefix) + "_R.txt", np.linalg.cholesky(Q).T)
    if (
        Path(str(prefix) + "_vectors.txt").exists()
        and "FLOATING_COMPLETE" in Path(str(prefix) + "_enum.log").read_text()
    ):
        result = subprocess.CompletedProcess([], 0)
    else:
        result = subprocess.run(
            [
                str(binary),
                str(n),
                str(minimum + 1e-7),
                "150000000",
                str(prefix) + "_R.txt",
                str(prefix) + "_vectors.txt",
                str(prefix) + "_enum.log",
            ],
            capture_output=True,
            text=True,
        )
    V = np.loadtxt(str(prefix) + "_vectors.txt", dtype=np.int64, ndmin=2)
    norm = np.einsum("vi,ij,vj->v", V, Q, V)
    assert np.all(norm == minimum), (name, np.unique(norm))
    M = V.T @ V
    L = np.linalg.cholesky(Q)
    e = np.linalg.eigvalsh(L.T @ M @ L)
    types = []
    for value in e:
        if not types or abs(value - types[-1][0]) > 1e-5:
            types.append([float(value), 1])
        else:
            types[-1][1] += 1
    row = {
        "name": name,
        "dimension": n,
        "source_minimum": minimum,
        "determinant": str(det),
        "source_center_density": (minimum / 4) ** (n / 2) / det**0.5,
        "returned_roots": len(V),
        "float_complete": result.returncode == 0,
        "moment_eigenvalues": types,
        "seconds": time.time() - start,
    }
    rows.append(row)
    print(json.dumps(row), flush=True)
    (P / "catalog_screen.json").write_text(json.dumps(rows, indent=2) + "\n")
