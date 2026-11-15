"""Exact finite-cardinality accounting for the three saved spherical proofs.

This reads accepted reports; it does not rerun their sign checkers, certify a
new density bound, or search for coefficients. Run from the project root.
"""
from fractions import Fraction as Q
import json
from pathlib import Path

ROOT = Path.cwd()
SOURCE = ROOT / "work/d4/radial-stratified-origin-bound/run-final/proof"
OUT = Path(__file__).with_suffix(".json")


def read(name):
    return json.loads((SOURCE / name).read_text())


def decimal(value):
    return float(value)  # Display only; every decision below uses Fraction.


def main():
    angular = [read(f"ANGULAR_{i}.json") for i in range(3)]
    local = read("LOCAL_ENUMERATION.json")
    acceptance = read("ACCEPTANCE.json")
    weight = Q(angular[0]["accepted_combination"]["positive_normalization"])
    weighted_B = Q(angular[0]["accepted_combination"]["weighted_accepted_B"])
    records = []
    for index, report in enumerate(angular):
        assert report["pair"]["status"] == report["triangle"]["status"] == "CERTIFIED"
        assert Q(report["accepted_combination"]["positive_normalization"]) == weight
        assert Q(report["accepted_combination"]["weighted_accepted_B"]) == weighted_B
        epsilon = Q(report["repair"]["scale"])
        eta = Q(report["repair"]["pair_shift"])
        H = Q(report["repair"]["diagonal"])
        delta = epsilon * H + eta
        M = report.get("maximum_code_count", report.get("fixed_count"))
        raw_B = Q(report["raw_B"])
        repaired_B = Q(report["normalized_B"] if index == 1 else report["B"])
        assert repaired_B == raw_B + (M - 1) * delta
        assert Q(report["repair"]["B_increment"]) == (M - 1) * delta
        assert epsilon > 0 and eta >= 0 and H > 0
        records.append(dict(index=index, M=M, epsilon=epsilon, eta=eta,
                            H=H, delta=delta, raw_B=raw_B, repaired_B=repaired_B))

    def applies(record, N, first_radius):
        if record["index"] == 1:
            return 0 <= N <= 26
        report = angular[record["index"]]
        return (N == report["fixed_count"] and
                first_radius <= report["condition_radius_index"] + 1)

    def gain(record, N):
        return weight * N * ((record["M"] - 1) * record["delta"]
                             - (N - 1) * record["eta"]
                             - record["epsilon"] * (N - 1) * (N - 2))

    profile = []
    for row in local["first_containing_radius_profile"]:
        N, first_radius = row["N"], row["first_containing_radius_index"]
        available = [r for r in records if applies(r, N, first_radius)]
        old_surpluses = [N * (weighted_B - weight * r["repaired_B"]) for r in available]
        new_surpluses = [s + gain(r, N) for r, s in zip(available, old_surpluses)]
        old_surplus = max([Q(0)] + old_surpluses)
        new_surplus = max([Q(0)] + old_surpluses + new_surpluses)
        assert old_surplus == Q(row["spherical_surplus"])
        assert new_surplus >= old_surplus
        new_floor = Q(row["exact_floor"]) + new_surplus - old_surplus
        profile.append(dict(N=N, first_containing_radius_index=first_radius,
                            old_floor=row["exact_floor"], new_floor=str(new_floor),
                            gain=str(new_surplus-old_surplus),
                            new_floor_decimal=decimal(new_floor),
                            minimizing_counts=row["minimizing_counts"]))
    assert len(profile) == 124
    assert {row["N"] for row in profile} == set(range(27))
    old_beta = min(Q(row["old_floor"]) for row in profile)
    new_beta = min(Q(row["new_floor"]) for row in profile)
    assert old_beta == Q(acceptance["new_origin_credit"])
    assert new_beta > old_beta
    raw = Q(acceptance["exact_source_objective"])
    target_b = 8 * (raw-Q(629, 1000))
    assert new_beta < target_b
    summary = dict(status="COMPUTED_EXACT_FINITE_COUNT_ACCOUNTING",
                   sign_checkers_rerun=False, density_bound_certified=False,
                   source_reports=[str(SOURCE.relative_to(ROOT)/f"ANGULAR_{i}.json") for i in range(3)],
                   weighted_gain_formula="W*N*((M-1)*delta-(N-1)*eta-epsilon*(N-1)*(N-2))",
                   input_profile="124 saved first-containing-radius minima; no new count-vector enumeration",
                   old_beta=str(old_beta), new_beta=str(new_beta),
                   new_beta_decimal=decimal(new_beta),
                   beta_gain=str(new_beta-old_beta), beta_gain_decimal=decimal(new_beta-old_beta),
                   target_b=str(target_b), remaining_energy_gap=str(target_b-new_beta),
                   remaining_energy_gap_decimal=decimal(target_b-new_beta),
                   derived_raw_objective=str(raw-new_beta/8),
                   derived_raw_objective_decimal=decimal(raw-new_beta/8),
                   minimizers=[r["minimizing_counts"] for r in profile if Q(r["new_floor"]) == new_beta],
                   certificate_gains=[dict(index=r["index"], N=N, gain=str(gain(r,N)),
                                           gain_decimal=decimal(gain(r,N)))
                                      for r in records for N in ([24,25,26] if r["index"]==1 else [r["M"]])],
                   profile=profile)
    OUT.write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k in
                      ["status", "new_beta_decimal", "beta_gain_decimal", "remaining_energy_gap_decimal",
                       "derived_raw_objective_decimal", "minimizers", "certificate_gains"]}, indent=2))


if __name__ == "__main__":
    main()
