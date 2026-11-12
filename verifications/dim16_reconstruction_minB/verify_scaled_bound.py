"""Exact transfer of the published 16D radial/count certificate to a new angular kernel.
Run verify_geometry.py first and separately replay original 16D verification.
"""
if not __debug__: raise RuntimeError('Assertions required')
from pathlib import Path
import json,hashlib
from fractions import Fraction as F
from flint import arb,fmpq,ctx
ctx.prec=256
D=Path(__file__).resolve().parent
BASE=D.parent/'dim16'
assert BASE.is_dir(), 'Required published base package missing'
def a(x):
 z=F(x);return arb(fmpq(z.numerator,z.denominator))
old=json.loads((BASE/'geometry'/'angular.json').read_text());new=json.loads((D/'angular.json').read_text())
rad=json.loads((BASE/'radial'/'candidate_exact.json').read_text())
assert (D/'fit.json').read_bytes()==(BASE/'geometry'/'fit.json').read_bytes()
check=json.loads((D/'RESULT_GEOMETRY.json').read_text());assert check['status']=='PASS'
for path,digest in check['input_hashes'].items():assert hashlib.sha256((D/path).read_bytes()).hexdigest()==digest
Bo=F(old['B_exact']);Bn=F(new['B_exact']);C=F(rad['C_exact']);Cn=C*Bn/Bo
assert Cn/Bn==C/Bo and 0<Cn<C
assert F(rad['normalization_exact'])==1
obj=F(rad['physical'][0])+Cn+F(rad['U_exact'])*F(rad['lambda_count_exact'])
density=arb.pi()**8/10321920*a(obj)
upper=F('0.024920883');assert density<a(upper)
out={'status':'PASS_CONDITIONAL_ON_ORIGINAL_RADIAL_COUNT_REPLAY','density_interval':str(density),'strict_upper':str(upper),'C_exact':str(Cn),'old_B_exact':str(Bo),'new_B_exact':str(Bn),'pair_function_unchanged':True,'original_dependency':'043b11949c629d2fe977db14f29ce174653c4068/verifications/dim16','scope':'New exact angular geometry and algebraic radial transfer. Original radial and count certificates must also be replayed.'}
(D/'TRANSFER_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
