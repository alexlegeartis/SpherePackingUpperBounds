"""Check the new finite-count consequence of the accepted spherical proofs.

No sign subdivision, optimization, or count census is run. Inherited certified
signs and grouped minima are dependencies, identified against their accepted
inputs. Run from the project root with the project Python, without -O.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
BASE = ROOT / "work/d4/radial-stratified-origin-bound"
SOURCE = BASE / "source"
PROOF = BASE / "run-final/proof"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    return json.loads(path.read_text(), object_pairs_hook=unique)


def stable(value):
    if isinstance(value, dict):
        return {k: stable(v) for k, v in value.items() if k != "seconds"}
    if isinstance(value, list):
        return [stable(v) for v in value]
    return value


def decimal(value):
    with localcontext() as context:
        context.prec = 40
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main():
    require(__debug__ and not sys.flags.optimize, "Python -O forbidden")
    acceptance = read(PROOF / "ACCEPTANCE.json")
    require(acceptance["status"] == "LOCAL_FULL_V7_AND_FRESH_STRATIFIED_ORIGIN_ACCEPTED"
            and acceptance["global_density_bound_proved"] is True,
            "missing accepted parent theorem")
    require(acceptance["count_vectors"] == 165503 and acceptance["radial_planes"] == 12
            and acceptance["angular_proofs_fresh"] == 3, "incomplete parent")

    # These are existing certificate-identity checks, not a new hash inventory.
    def bind(path):
        key = str(path.relative_to(ROOT))
        require(key in acceptance["source_bindings"], "unbound inherited input")
        require(hashlib.sha256(path.read_bytes()).hexdigest()
                == acceptance["source_bindings"][key], "changed inherited input: " + key)
        return path

    registry = read(bind(SOURCE / "target625/origin_completion/FINAL_ORIGIN_REGISTRY.json"))
    algebra_path = bind(SOURCE / "target625/geometry/fixed_cardinality_spherical_algebra.py")
    spec = importlib.util.spec_from_file_location("finite_count_algebra", algebra_path)
    algebra = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(algebra)
    for name in ("certify_fixedN_common_spherical.py", "certify_accepted_combined_boundedN_credit.py",
                 "certify_accepted_stratified_fixedN_credit.py"):
        bind(SOURCE / "target625/geometry" / name)
    origin = read(bind(SOURCE / "current/candidates/round9_six_graphs_degree18/origin.json"))
    rays = origin["appended_spherical_local_energy_rays"]
    W = sum(Q(ray["weight"]) for ray in rays)
    Bbar = sum(Q(ray["weight"]) * Q(ray["B"]) for ray in rays)
    gbar = [sum(Q(ray["weight"]) * Q(ray["scalar_coefficients"][j]) for ray in rays)
            for j in range(13)]
    records = []
    for index, entry in enumerate(registry["common_g_certificates"]):
        report = read(PROOF / f"ANGULAR_{index}.json")
        prior = read(bind(SOURCE / entry["prior_verification"]))
        require(stable(report) == stable(prior), "report differs from accepted source")
        require(report["pair"]["status"] == report["triangle"]["status"] == "CERTIFIED",
                "missing inherited sign proof")
        data = read(bind(SOURCE / entry["certificate"]))
        require(Q(data["feature_scale"]) > 0, "invalid Gram scale")
        matrices = algebra.gram_blocks(data)  # Checks integer rays and fixed-N restriction.
        raw_B = Q(sum(data["scalar"]), data["scalar_scale"]) + Q(
            int(matrices[0][0, 0]), data["feature_scale"] ** 2)
        require(raw_B == Q(report["raw_B"]), "raw origin identity")
        common = report["accepted_combination"]
        require(Q(common["positive_normalization"]) == W
                and Q(common["weighted_accepted_B"]) == Bbar
                and list(map(Q, common["weighted_accepted_g"])) == gbar
                and list(map(Q, report["common_scalar_coefficients"])) == [x/W for x in gbar],
                "physical spherical energy mismatch")
        epsilon, eta = map(Q, (report["repair"]["scale"], report["repair"]["pair_shift"]))
        H = Q(report["repair"]["diagonal"])
        delta = epsilon*H + eta
        M = report.get("maximum_code_count", report.get("fixed_count"))
        old_B = Q(report["normalized_B"] if index == 1 else report["B"])
        require(epsilon > 0 and eta >= 0 and old_B == raw_B+(M-1)*delta,
                "repair accounting mismatch")
        if index == 1:
            require(M == 26 and not data.get("fixed_count") and Q(data["s"]) == Q(21,40),
                    "bounded scope mismatch")
        else:
            require(data["fixed_count"] == M == (24 if index == 0 else 25), "fixed count")
            R = Q(report["condition_radial_radius"])
            require(R*R < 2 and 1-Q(1)/(2*R*R) <= Q(data["s"]), "radius scope")
        records.append(dict(index=index, M=M, epsilon=epsilon, eta=eta, raw_B=raw_B,
                            old_B=old_B, first=report.get("condition_radius_index",4)+1))

    local = read(PROOF / "LOCAL_ENUMERATION.json")
    rows = local["first_containing_radius_profile"]
    require(len(rows) == 124 and len({(r["N"],r["first_containing_radius_index"]) for r in rows}) == 124,
            "incomplete or duplicate grouped profile")
    require({r["N"] for r in rows} == set(range(27)), "missing cardinality")
    updated = []
    for row in rows:
        N, first = row["N"], row["first_containing_radius_index"]
        available = [r for r in records if r["index"] == 1 or (N == r["M"] and first <= r["first"])]
        old = max([Q(0)] + [N*(Bbar-W*r["old_B"]) for r in available])
        direct = [N*Bbar-W*(N*r["raw_B"]+r["eta"]*N*(N-1)
                          +r["epsilon"]*N*(N-1)*(N-2)) for r in available]
        new = max([old]+direct)
        require(old == Q(row["spherical_surplus"]), "saved profile inconsistent")
        floor = Q(row["exact_floor"])+new-old
        updated.append(dict(N=N, first_radius=first, old_floor=row["exact_floor"],
                            new_floor=str(floor), gain=str(new-old),
                            minimizing_counts=row["minimizing_counts"]))
    old_beta = min(Q(row["old_floor"]) for row in updated)
    beta = min(Q(row["new_floor"]) for row in updated)
    require(old_beta == Q(acceptance["new_origin_credit"]) < beta, "no improvement")
    candidate = read(HERE / "multiplicity/finite_count_gain.json")
    require(beta == Q(candidate["new_beta"]), "independent arithmetic mismatch")
    raw = Q(acceptance["exact_source_objective"])
    require(raw-old_beta/8 == Q(acceptance["objective"]), "parent normalization")
    graph_saving = Q(3,8000000)  # Previously proved capacity 1358 -> 1208.
    objective = raw-graph_saving-beta/8
    rounded = Q((objective*10**12).numerator//(objective*10**12).denominator+1,10**12)
    old_ceiling = Q(125907730731,200000000000)
    b = 8*(raw-Q(629,1000))
    require(objective < rounded < old_ceiling and beta < b, "target or rounding mismatch")
    report = dict(status="EXACT_FINITE_COUNT_COMPOSITION_ACCEPTED", dimension=4,
                  packing_family="all packings, including arbitrary upper density",
                  evidence="PROVED finite-count consequence of inherited certified signs; COMPUTED EXACT objective",
                  inherited_sign_checks_rerun=False, new_count_census=False,
                  source= str(PROOF.relative_to(ROOT)), grouped_minima_checked=124,
                  original_count_vectors=165503, exact_raw_objective=str(raw),
                  old_beta=str(old_beta), new_beta=str(beta), new_beta_decimal=decimal(beta),
                  beta_gain=str(beta-old_beta), graph_saving=str(graph_saving),
                  exact_objective=str(objective), objective_decimal=decimal(objective),
                  strict_upper_bound=str(rounded), strict_upper_bound_decimal=decimal(rounded),
                  previous_strict_upper_bound=str(old_ceiling),
                  displayed_improvement=str(old_ceiling-rounded),
                  target_b=str(b), remaining_energy_gap=str(b-beta),
                  remaining_energy_gap_decimal=decimal(b-beta),
                  seven_classes=[r for r in updated if Q(r["new_floor"]) < b],
                  minimizers=[r["minimizing_counts"] for r in updated if Q(r["new_floor"]) == beta],
                  profile=updated)
    require(len(report["seven_classes"]) == 7, "changed poor-class count")
    output = HERE / "FINITE_COUNT_ACCEPTANCE.json"
    temp = output.with_suffix(".tmp")
    temp.write_text(json.dumps(report,indent=2)+"\n")
    temp.replace(output)
    print(json.dumps({k:report[k] for k in ("status","new_beta_decimal","objective_decimal",
                      "strict_upper_bound_decimal","remaining_energy_gap_decimal")},indent=2))


if __name__ == "__main__":
    main()
