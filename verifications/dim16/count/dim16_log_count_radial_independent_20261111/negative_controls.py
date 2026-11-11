"""Explicit adversarial controls for the independent radial gate."""
if not __debug__:raise RuntimeError('Assertions required')
import verify as V
import json,subprocess,sys
from fractions import Fraction as F
from flint import arb,arb_poly,fmpq as Q
results={}
def rejected(name,fun):
 try:fun()
 except (AssertionError,RuntimeError):results[name]={'rejected':True}
 else:raise RuntimeError('Corruption accepted: '+name)

def validate_cover(rows):
 rows=sorted(rows,key=lambda z:F(z['lo']))
 assert F(rows[0]['lo'])==1 and F(rows[-1]['hi'])==V.end
 for i,z in enumerate(rows):
  lo,hi=F(z['lo']),F(z['hi']);assert lo<hi
  if i:assert F(rows[i-1]['hi'])==lo
  assert not(lo<V.threshold<hi)
  assert bool(z['count_active'])==((lo+hi)/2<V.threshold)
  assert all(not(lo<k<hi)for k in V.VK+V.LK)
 assert any(F(z['lo'])==V.threshold and not z['count_active']for z in rows)
validate_cover(V.saved)
rejected('missing_pair_interval',lambda:validate_cover(V.saved[1:]))
wrong=[dict(z,count_active=True)if F(z['lo'])==V.threshold else z for z in V.saved]
rejected('count_allowance_wrongly_retained_on_right',lambda:validate_cover(wrong))
# Wrong Fourier eigenvalue in the epsilon repair cannot satisfy the direct transform.
wrongH=(V.Hs/V.N+V.eps*V.E)/(1-V.eps)
assert wrongH!=V.Hdirect and wrongH(0)!=1
results['wrong_fourier_repair_sign']={'rejected':True}
# The original numerical witness really violated a full pair gate. This is not
# just a metadata rejection: interval arithmetic proves a positive point.
xx=V.aq(F(6858771052618467,562949953421312));v=xx/V.TP
assert v>V.aq(F(7,4)) and v<V.aq(F(2))
pa=V.parts(F(39,20));phi=V.phival(v,pa);h=(-phi).exp()
oldP=arb_poly([arb(z)for z in V.Ps.coeffs()]);oldC=arb(Q(V.source['C']))
residual=oldP(xx)*(-xx).exp()-oldC/V.aq(V.B)*h*(2-h)
assert residual>0
results['unrepaired_source_rejected_by_positive_pair_point']={'rejected':True,'X_exact':'6858771052618467/562949953421312','physical_residual_interval':str(residual)}
# A too-small claimed rounded result must fail the rational upper gate.
assert V.dlow>F('0.02492098')
results['overstated_decimal_bound']={'rejected':True}
r=subprocess.run([sys.executable,'-O',str(V.D/'verify.py')],capture_output=True,text=True)
assert r.returncode!=0 and 'Run without optimization'in r.stderr
results['python_optimization_disabled']={'rejected':True}
(V.D/'NEGATIVE_CONTROLS.json').write_text(json.dumps({'status':'PASS','controls':results},indent=2)+'\n')
print('NEGATIVE CONTROLS PASS',results)
