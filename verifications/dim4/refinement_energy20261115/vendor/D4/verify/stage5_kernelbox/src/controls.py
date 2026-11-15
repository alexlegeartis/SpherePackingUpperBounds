"""Lemma 1123 (kernelbox-package): negative and malformed-input controls.

Each control feeds check_kernelbox.py (or kbox_run.py) a wrong input and passes
only if the tool fails closed (nonzero exit, status FAILED) for the expected reason.
Inputs and outputs go to paths.OUT/controls/. Run from any folder:
  nice -n 10 python3 src/controls.py
"""
from pathlib import Path
import copy, json, subprocess, sys, time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import paths
C = paths.OUT / "controls"
C.mkdir(parents=True, exist_ok=True)
PY = sys.executable
j1 = json.loads((paths.KBOX / "kbox_j1.json").read_text())
DIP = 27  # root t-chunk [-0.8656, -0.8606] containing the dip


def write(name, obj=None, text=None):
    p = C / name
    p.write_text(text if text is not None else json.dumps(obj))
    return p


inputs = {}
d = copy.deepcopy(j1); d["summary"]["lam"] = "4/1000000"
inputs["neg_lam4"] = (write("in_neg_lam4.json", d), ["--lam1", "4/1000000"], "negative leaf")
d = copy.deepcopy(j1); d["roots"][DIP]["tree"] = "L"
inputs["pruned_dip_tree"] = (write("in_pruned_dip.json", d), [], "negative leaf")
d = copy.deepcopy(j1); del d["roots"][5]
inputs["missing_root"] = (write("in_missing_root.json", d), [], "roots do not match")
d = copy.deepcopy(j1); d["roots"][DIP]["tree"] = [7, "L", "L"]
inputs["bad_split_dim"] = (write("in_bad_dim.json", d), [], "bad node")
inputs["truncated_json"] = (write("in_truncated.json", text=(paths.KBOX / "kbox_j1.json").read_text()[:5000]), [], "JSONDecodeError")
d = copy.deepcopy(j1); d["summary"]["status"] = "FAILED"
inputs["failed_status"] = (write("in_failed_status.json", d), [], "AssertionError")

results = {}
for name, (path, extra, expect) in inputs.items():
    out = C / f"out_{name}.json"
    t0 = time.time()
    r = subprocess.run([PY, str(HERE / "check_kernelbox.py"), "--workers", "4", "--j1-file", str(path),
                        "--out", str(out)] + extra, capture_output=True, text=True)
    rec = json.loads(out.read_text())
    ok = r.returncode != 0 and rec["status"] == "FAILED" and expect in rec.get("error", "")
    results[name] = dict(exit=r.returncode, status=rec["status"], error=rec.get("error", "")[:300],
                         expected=expect, control_pass=ok, seconds=round(time.time() - t0, 1))
    print(name, "PASS" if ok else "UNEXPECTED", r.returncode, rec.get("error", "")[:120], flush=True)

# independent negative control: the search itself must fail at lam = 4/10^6
t0 = time.time()
r = subprocess.run([PY, str(HERE / "kbox_run.py"), "--j", "1", "--lam", "4/1000000", "--workers", "4",
                    "--out", str(C / "kbox_j1_lam4.json")], capture_output=True, text=True)
summ = json.loads((C / "kbox_j1_lam4.json").read_text())["summary"]
ok = r.returncode == 1 and summ["status"] == "FAILED" and len(summ["failures"]) > 0
results["search_lam4"] = dict(exit=r.returncode, status=summ["status"], nfail=len(summ["failures"]),
                              first_fail=summ["failures"][:1], control_pass=ok, seconds=round(time.time() - t0, 1))
print("search_lam4", "PASS" if ok else "UNEXPECTED", r.returncode, len(summ["failures"]), flush=True)

allok = all(v["control_pass"] for v in results.values())
(C / "controls.json").write_text(json.dumps(dict(status="ALL_CONTROLS_FAIL_CLOSED" if allok else "CONTROL_PROBLEM",
                                                 controls=results), indent=1) + "\n")
print("ALL_CONTROLS_FAIL_CLOSED" if allok else "CONTROL_PROBLEM")
sys.exit(0 if allok else 1)
