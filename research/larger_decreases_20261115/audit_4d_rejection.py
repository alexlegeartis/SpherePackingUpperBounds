"""Find and exactly verify a cyclic-sign violation in the saved 4D proposal."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,sys
import numpy as np
from search_spherical_energy import Model

if not __debug__:raise RuntimeError('assertions must be enabled')
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
proposal=np.load(HERE/'dim4_research/joint_energy/proposal.npz')
n,d,L,D=[int(proposal[k]) for k in ['N','radial_degree','angular_degree','p_degree']]
s=Q(float(proposal['s']));model=Model(n,float(s),d,L,D)
values=proposal['values']
grid=[Q(k,16) for k in range(-16,9)]
points=[z for z in itertools.combinations_with_replacement(grid,3)
        if 1+2*z[0]*z[1]*z[2]-sum(x*x for x in z)>=0]
numeric=np.array([[float(x) for x in z] for z in points])
scores=model.triangle(numeric)@values
order=np.argsort(scores)[::-1]
matrices=[[[Q(float(x)) for x in row] for row in matrix] for matrix in proposal['matrices']]
free=list(map(lambda x:Q(float(x)),proposal['free']))
def basis(t):
    x=2*(t+1)/(1+s)-1
    out=[Q(1),x]
    for k in range(2,d+1):out.append(2*x*out[-1]-out[-2])
    return out
def kernel(a,b,t):
    ba,bb=basis(a),basis(b)
    value=free[1]+sum(free[i+1]*(ba[i]+bb[i]) for i in range(1,d+1))
    value+=sum(ba[i]*matrices[0][i][j]*bb[j] for i in range(1,d+1) for j in range(1,d+1))
    x=t-a*b;w=(1-a*a)*(1-b*b);h0,h1=Q(1),x
    for ell in range(1,L+1):
        value+=h1*sum(ba[i]*matrices[ell][i][j]*bb[j] for i in range(d+1) for j in range(d+1))
        h0,h1=h1,((2*ell+1)*x*h1-ell*w*h0)/(ell+1)
    return value/Q(10**6)
for index in order[:40]:
    a,b,t=points[index]
    value=kernel(a,b,t)+kernel(a,t,b)+kernel(b,t,a)
    if value>0:
        report=dict(status='REJECTED_EXACT_POSITIVE_CYCLIC_VALUE',certified_density_bound=False,
                    point=list(map(str,[a,b,t])),gram_determinant=str(1+2*a*b*t-a*a-b*b-t*t),
                    cyclic_value_exact=str(value),cyclic_value_decimal=float(value),
                    interpretation='The saved binary64 matrices and free coefficients are treated as exact rational numbers; the cyclic sign fails on a valid spherical triple.')
        (HERE/'dim4_research/EXACT_REJECTION.json').write_text(json.dumps(report,indent=2)+'\n')
        print('REJECTED: exact positive cyclic value',float(value),'at',list(map(str,[a,b,t])))
        break
else:raise RuntimeError('No exact violating point found on the diagnostic grid')
