"""Replay only the recovered fixed graph certificates and local helper audit.

Run from the project root with the project Python. No optimization, global
pair/cyclic replay, or writes to the retained source occur here. Outputs go to
a fresh directory so historical review reports are preserved.
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
# Stage 1 export: the record and output directory come from the command line
# (--record, --output-dir); subprocesses use this interpreter and run in the
# record. The independent moderate checker sits next to this file.
HERE = Path(__file__).resolve().parent
RECORD = None
ROUTE = None
OUTPUT = None
PYTHON = Path(sys.executable)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save(name, data):
    target = OUTPUT / name
    require(not target.exists(), f"Refusing to overwrite an existing report: {target}")
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(target)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def command_report(checker, certificate, name, status, capacity):
    args = [str(PYTHON), "-B", str(checker), str(certificate)]
    run = subprocess.run(args, cwd=RECORD, capture_output=True, text=True,
                         timeout=60, check=True)
    require(not run.stderr, f"Unexpected stderr: {run.stderr}")
    report = json.loads(run.stdout)
    require(report["status"] == status, "Unexpected verifier status")
    require(report["even_ordered_edge_bound"] == capacity, "Capacity mismatch")
    require(report["global_density_bound_proved"] is False, "Wrong graph scope")
    save(name, report)
    return dict(command=args, returncode=run.returncode, output=name,
                expected_status=status, even_ordered_edge_bound=capacity)


def controls(verifier, certificate):
    original = json.loads(certificate.read_text())
    tests = []

    def rejected(name, mutation):
        candidate = copy.deepcopy(original)
        mutation(candidate)
        try:
            verifier.verify(candidate)
        except (AssertionError, ValueError, KeyError, TypeError, ZeroDivisionError) as error:
            tests.append(dict(name=name, result="REJECTED", error_type=type(error).__name__,
                              diagnostic=str(error)))
        else:
            raise ValueError(f"Malformed certificate was accepted: {name}")

    rejected("missing_power_coefficients", lambda d: d.pop("power_coefficients"))
    rejected("incomplete_harmonic_vector", lambda d: d["harmonic_coefficients"].pop())
    rejected("negative_harmonic_coefficient", lambda d: d["harmonic_coefficients"].__setitem__(1, "-1"))
    rejected("wrong_power_identity", lambda d: d["power_coefficients"].__setitem__(0, "0"))
    rejected("wrong_even_capacity_1878", lambda d: d.__setitem__("claimed_even_ordered_edge_bound", 1878))
    rejected("wrong_ambient_dimension", lambda d: d.__setitem__("ambient_dimension", 5))
    rejected("wrong_separation", lambda d: d.__setitem__("minimum_pair_distance_squared", "2"))
    rejected("nonintegral_harmonic_degree", lambda d: d.__setitem__("harmonic_degree", 16.5))
    rejected("malformed_rational", lambda d: d["harmonic_coefficients"].__setitem__(0, "not-a-rational"))
    rejected("noncanonical_rational", lambda d: d["harmonic_coefficients"].__setitem__(10, "0/2"))
    rejected("numeric_float_coefficient", lambda d: d["harmonic_coefficients"].__setitem__(10, 0.0))
    return tests


def helper_replay():
    # The archived audit contains the original absolute scratch paths. Rebind
    # those two path fields in a staged copy; preserve every source byte and
    # every original identity hash. No archived writer is invoked.
    raw = json.loads((RECORD / "geometry/SIX_GRAPHS_SMALL_CAP_BOUNDS_HELPER_CONSTRUCTION.json").read_text())
    historical_paths = {key: raw[key] for key in ("original_path", "new_path")}
    for key in historical_paths:
        raw[key] = str(RECORD / raw[key].split("/93cb6ce085de/", 1)[1])
    save("helper_construction_local_paths.json", raw)
    module = load("_intake_graph_archived_helper_audit",
                  RECORD / "geometry/audit_six_graphs_small_cap_bounds_helper.py")
    module.ROOT = RECORD
    module.HERE = RECORD / "geometry"
    module.ROUTE = ROUTE
    module.CONSTRUCTION = OUTPUT / "helper_construction_local_paths.json"
    report = module.audit()
    require(report["status"] == "EXACT_SIXTH_GRAPH_PROOF_AND_HELPER_BOOKKEEPING_AUDIT_PASS",
            "Unexpected helper audit status")
    require(report["singleton_controls"] == 1406, "Unexpected singleton count")
    require(sum(row["exact_direct_point_controls"] for row in report["box_controls"]) == 270,
            "Unexpected box point control count")
    require(report["default_zero_exactly_recovers_parent"] is True, "Parent recovery fails")
    require(report["full_branches_exactly_unchanged"] is True, "Full branches changed")
    require(report["global_density_bound_proved"] is False, "Wrong audit scope")
    report["intake_replay_path_rebinding"] = historical_paths
    save("helper_bookkeeping_replay.json", report)
    return report


def main():
    global RECORD, ROUTE, OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, required=True,
                        help="Working copy of the record.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Fresh report directory; no existing output file is overwritten.")
    options = parser.parse_args()
    require(__debug__, "Assertions must remain enabled")
    RECORD = options.record.resolve()
    ROUTE = RECORD / "current/six_exterior_graph_route"
    OUTPUT = options.output_dir.resolve()
    require(RECORD.is_dir(), "Record is not present")
    require(OUTPUT != RECORD and RECORD not in OUTPUT.parents, "Outputs cannot modify the record")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    require(not any(OUTPUT.iterdir()), "Choose a fresh output directory; prior reports are preserved")
    started = time.time()
    commands = []
    for checker, certificate, name, status, capacity in [
        ("verify_moderate_exterior_graph.py", "moderate_exterior_graph_certificate.json", "moderate_bernstein.json",
         "EXACT_SINGLE_COLOR_EXTERIOR_GRAPH_CERTIFICATE_PASS", 1880),
        ("verify_outer_exterior_graph.py", "outer_exterior_graph_certificate.json", "near_antipodal_120.json",
         "EXACT_NEAR_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS", 120),
        ("verify_antipodal_wide_exterior_graph.py", "antipodal_wide_exterior_graph_certificate.json", "broad_antipodal_1358.json",
         "EXACT_PROJECTED_S2_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS", 1358),
        ("verify_narrow_exterior_graph.py", "narrow_exterior_graph_certificate.json", "narrow_antipodal_604.json",
         "EXACT_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS", 604),
    ]:
        commands.append(command_report(ROUTE / checker, ROUTE / certificate,
                                       name, status, capacity))

    independent_checker = HERE / "verify_moderate_independent.py"
    args = [str(PYTHON), "-B", str(independent_checker),
            str(ROUTE / "moderate_exterior_graph_certificate.json")]
    result = subprocess.run(args, cwd=RECORD, capture_output=True, text=True,
                            timeout=60, check=True)
    require(not result.stderr, f"Unexpected independent stderr: {result.stderr}")
    independent = json.loads(result.stdout)
    require(independent["status"] == "INDEPENDENT_EXACT_MODERATE_GRAPH_CAPACITY_1880_PASS",
            "Unexpected independent checker status")
    require(independent["even_ordered_edge_capacity"] == 1880, "Independent capacity mismatch")
    require([row["interior_root_count"] for row in independent["whole_interval_sturm_proofs"]] == [0, 0],
            "Independent root count mismatch")
    save("moderate_independent_sturm.json", independent)
    commands.append(dict(command=args, returncode=result.returncode,
                         output="moderate_independent_sturm.json",
                         expected_status=independent["status"]))

    sys.path.insert(0, str(ROUTE))
    moderate = load("_intake_graph_moderate_checker", ROUTE / "verify_moderate_exterior_graph.py")
    negative = controls(moderate, ROUTE / "moderate_exterior_graph_certificate.json")
    save("moderate_negative_controls.json", dict(status="ALL_EXPECTED_REJECTIONS", controls=negative))
    helper = helper_replay()
    manifest_path = RECORD / "current/candidates/round9_six_graphs_degree18/graph_energies.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest["schema"] == "INDEPENDENT_EXTERIOR_GRAPH_ENERGIES_V7", "Wrong graph manifest schema")
    weights = tuple(map(Q, manifest["weights"]))
    near_weight, narrow_weight, moderate_weight = (
        Q(manifest[key]) for key in ("outer_graph_weight", "narrow_graph_weight", "moderate_graph_weight"))
    require(weights == (Q(0), Q(0), Q(1,50000000)), "Unexpected V7 A/B/E9 weights")
    require((near_weight, narrow_weight, moderate_weight) ==
            (Q(1,5000000), Q(1,5000000), Q(1,10**8)), "Unexpected active graph weights")
    # These hashes are the producer's explicit graph source-binding contract,
    # not a new inventory of unrelated files.
    source_bindings = {}
    for relative, expected in manifest["source_sha256"].items():
        require(Path(relative).name == relative, "Graph binding escaped route")
        actual = hashlib.sha256((ROUTE / relative).read_bytes()).hexdigest()
        require(actual == expected, f"Graph manifest source mismatch: {relative}")
        source_bindings[relative] = actual
    charges = [Q(1358)*weights[2]/8, Q(120)*near_weight/8,
               Q(604)*narrow_weight/8, Q(1880)*moderate_weight/8]
    require(sum(charges) == Q(4769,200000000), "Graph charge mismatch")
    save("graph_manifest_binding.json", dict(status="EXACT_GRAPH_SOURCE_MANIFEST_AND_WEIGHTS_PASS",
         manifest=str(manifest_path.relative_to(RECORD)), source_sha256=source_bindings,
         weights=list(map(str,weights)), outer_graph_weight=str(near_weight),
         narrow_graph_weight=str(narrow_weight), moderate_graph_weight=str(moderate_weight),
         combined_density_charge=str(sum(charges)), global_density_bound_proved=False))
    summary = dict(status="GRAPH_LOCAL_ACCEPTANCE_PASS", record=str(RECORD),
                   output_dir=str(OUTPUT), commands=commands,
                   malformed_input_rejections=len(negative), singleton_controls=helper["singleton_controls"],
                   rational_box_point_controls=270,
                   combined_graph_density_charge=str(sum(charges)),
                   graph_density_charges=list(map(str,charges)),
                   global_density_bound_proved=False, elapsed_seconds=time.time()-started)
    save("ACCEPTANCE.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
