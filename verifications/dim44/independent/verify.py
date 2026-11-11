"""Portable replay of the frozen 44D audit and its independent checks."""
import hashlib,json,subprocess,sys
from pathlib import Path
if not __debug__: raise RuntimeError("Assertions must be enabled; do not use -O")
D=Path(__file__).resolve().parent
for name,expected in json.loads((D/"INPUT_SHA256.json").read_text()).items():
    actual=hashlib.sha256((D/name).read_bytes()).hexdigest()
    assert actual==expected,(name,actual,expected)
for script,log in [("frozen_package/verify.py","complete_replay.log"),("independent_radial.py","independent_radial.log")]:
    result=subprocess.run([sys.executable,str(D/script)],cwd=D,capture_output=True,text=True)
    (D/log).write_text(result.stdout+result.stderr)
    print(result.stdout,flush=True)
    if result.stderr: print(result.stderr,file=sys.stderr,flush=True)
    assert result.returncode==0,(script,result.returncode)
primary=json.loads((D/"frozen_package/full_replay_result.json").read_text())
independent=json.loads((D/"independent_radial.json").read_text())
assert primary["pass"] and primary["strict_upper"]==independent["strict_upper"]=="0.000001842868599"
for key in ["independent_eigenbasis_fourier","independent_orthogonality_PSD","origin_rank_one","mean_value_pair_pass","exact_rational_rounding"]:
    assert independent[key] is True,key
assert independent["source_sha256"]=="4ca48c456c914770630100c4f3865967a5a6e039c7a29d78ebadaca2b77c5c55"
out={"pass":True,"dimension":44,"strict_upper":"0.000001842868599","primary":primary,"independent":independent}
(D/"audit_result.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2),flush=True)
