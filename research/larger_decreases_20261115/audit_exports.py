"""Reconstruct each rational witness from its retained discovery coordinates."""

import argparse, hashlib, json, sys
from fractions import Fraction
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P

sys.path.insert(
    0, str(Path(__file__).resolve().parent)
)
from common import ROOT, SOURCE, source

HERE = Path(__file__).resolve().parent
if not __debug__:
    raise RuntimeError("Assertions must remain enabled")
sys.set_int_max_str_digits(0)
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", type=int, nargs="+")
args = ap.parse_args()
results = []
for n in args.dimensions:
    D = ROOT / f"verifications/dim{n}/refinement_larger20261115"
    r = json.loads((D / "certificate.json").read_text())
    candidate = ROOT / r["research_source"]
    candidate = candidate.with_name("candidate.json")
    c = json.loads(candidate.read_text())
    d = c["degree"]
    rawp, rawh, basis = source(n, d)
    coefficients = [Q(str(Fraction(float(z)))) for z in c["x"]]
    p = coefficients[0] * rawp + sum(
        (coefficients[k + 1] * basis[k] for k in range(d + 1)), P([])
    )
    h = coefficients[0] * rawh + sum(
        ((-1) ** k * coefficients[k + 1] * basis[k] for k in range(d + 1)), P([])
    )
    assert p == P(list(map(Q, r["physical_polynomial"])))
    assert h == P(list(map(Q, r["fourier_polynomial"])))
    assert coefficients[-1] * Q(r["C_guard_exact"]) == Q(r["C_exact"])
    assert d == r["correction_degree"] == 95 and p.degree() == h.degree() == 111
    for key in [
        "dimension",
        "s_exact",
        "a_exact",
        "b_exact",
        "q_exact",
        "B_exact",
        "angular_file",
    ]:
        assert c[key] == r[key], (n, key)
    assert (p(0) + Q(r["C_exact"])) / h(0) == Q(r["objective_exact"])
    angular = json.loads((ROOT / c["angular_file"]).read_text())
    frozen = json.loads((D / "angular.json").read_text())
    for key in ["gegenbauer_coefficients", "polynomial_coefficients"]:
        while Fraction(angular[key][-1]) == 0:
            angular[key].pop()
        assert angular[key] == frozen[key]
    assert angular["B_exact"] == frozen["B_exact"] == r["B_exact"]
    result = {
        "dimension": n,
        "pass": True,
        "candidate": str(candidate.relative_to(ROOT)),
        "candidate_sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
        "binary64_coordinates_exported_exactly": True,
        "both_polynomials_reconstructed_exactly": True,
        "origin_guard_and_objective_exact": True,
        "angular_data_exact_up_to_zero_padding": True,
        "correction_degree": d,
        "physical_degree": p.degree(),
    }
    results.append(result)
    print("EXPORT PASS", n, flush=True)
out = {
    "pass": True,
    "source_48D_seed": str(SOURCE.relative_to(ROOT)),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "results": results,
    "scope": "Provenance/export check; complete angular and radial verifications are separately required",
}
(HERE / "EXACT_EXPORT_AUDIT.json").write_text(json.dumps(out, indent=2) + "\n")
