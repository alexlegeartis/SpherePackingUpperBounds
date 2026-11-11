"""Build/replay exact rational physical-tail budgets for long-edge lens cuts.

Construction picks 85% of a numerical minimum, then proves every whole bin
with exact source coefficients and 256-bit Arb Taylor bounds. The independent
replay mode ignores construction minima and recomputes every sign proof.
No kernel or density bound is certified by this program.
"""
if not __debug__:
    raise RuntimeError('Do not disable assertions.')

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from flint import fmpq as Q, fmpq_poly as Poly, arb, ctx

ctx.prec = 256
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCE = HERE.parent / 'density_low_structural_20261111'
RADIUS = F(27,25)
CAP_THRESHOLDS = [(F(41,25),83),(F(17,10),8),(F(7,4),5),
                  (F(9,5),3),(F(46,25),2),(F(48,25),1)]
SOURCE_NAMES = {5:'patch_dim5_v2_exact.json',6:'patch5_exact.json',7:'patch_dim7_exact.json'}


def ev(p, x):
    value = arb(0)
    for c in reversed(p.coeffs()):
        value = value*x + arb(c)
    return value


def mag(p, x, h):
    co = [arb(0)]*(p.degree()+1)
    for c in reversed(p.coeffs()):
        for j in range(p.degree(),0,-1):
            co[j] = co[j]*x + co[j-1]
        co[0] = co[0]*x + arb(c)
    value = arb(0)
    for c in reversed(co):
        value = value*h + abs(c)
    return value


def geometry(n=None,mode='centroid',spherical_certificates=None):
    assert mode in ['centroid','midpoint_rankin','transverse_rankin','spherical_transverse']
    records = []
    if mode=='spherical_transverse':
        from spherical_lens_caps import projection,exact_certificate
        if spherical_certificates is None:
            all_codes=json.loads((HERE/'spherical_lens_caps.json').read_text())['certificates']
            spherical_certificates=[row for row in all_codes if row['packing_dimension']==n
                                    and row['projection']=='transverse' and F(row['base_lower'])<F(41,25)]
        assert [F(row['base_lower'])for row in spherical_certificates]==[F(3,2),F(31,20),F(8,5)]
        for row in spherical_certificates:
            expected=projection(n,F(row['base_lower']),'transverse')
            assert all(row[k]==v for k,v in expected.items())
            certificate=exact_certificate(row['ambient_dimension'],row['gegenbauer_coefficients'],Q(row['max_inner_product']))
            assert all(row[k]==v for k,v in certificate.items())
            records.append(dict(lower_base=row['base_lower'],cap=row['integer_bound'],spherical_certificate=row))
    for d,B in CAP_THRESHOLDS:
        rho2 = RADIUS**2-d**2/4
        bound = 1/(1-2*rho2)
        assert 0 <= rho2 < F(1,2) and B <= bound < B+1
        row=dict(lower_base=str(d),cap=B,midpoint_radius_squared=str(rho2),
                 real_count_upper=str(bound))
        if mode!='centroid':
            assert n>=2
            dot=rho2-F(1,2)
            dim=n
            if mode in ['transverse_rankin','spherical_transverse']:
                dot+=(RADIUS**2-1)**2/(2*d*d)
                dim=n-1
            assert dot<0
            row.update(centroid_cap=B,cap=min(B,dim+1),
                       obtuse_vector_dimension=dim,obtuse_vector_dot_upper=str(dot),
                       obtuse_cap=dim+1)
        records.append(row)
    return records


def cap_at(d,n=None,mode='centroid'):
    return cap_from_geometry(d,geometry(n,mode))


def cap_from_geometry(d,records):
    assert F(records[0]['lower_base']) <= d <= 2*RADIUS
    return [row['cap']for row in records if F(row['lower_base'])<=d][-1]


class PhysicalFunction:
    def __init__(self,path=None,payload=None):
        self.path = path
        if path is not None:
            self.sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
            source = json.loads(path.read_text())[0]
            self.n = source['dimension']
            self.source = source
            coeff = list(map(lambda c:str(F(str(c))),source['coefficients']))
            self.payload = dict(dimension=self.n,basis='L_k^(n/2-1)(2X)/L_k^(n/2-1)(0)',
                                variable='X=pi*r^2',physical_factor='exp(-X)',
                                additional_physical_terms='none',degree=len(coeff)-1,
                                exact_rational_coefficients=coeff)
        else:
            self.payload = payload
            self.n = payload['dimension']
            self.sha256 = None
            coeff = payload['exact_rational_coefficients']
            assert payload['basis']=='L_k^(n/2-1)(2X)/L_k^(n/2-1)(0)'
            assert payload['variable']=='X=pi*r^2' and payload['physical_factor']=='exp(-X)'
            assert payload['additional_physical_terms']=='none' and payload['degree']==len(coeff)-1
            assert coeff==[str(F(c))for c in coeff]
        self.physical_sha256 = hashlib.sha256(json.dumps(self.payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        X = Poly([0,1]); m = Q(self.n,2)
        ls = [Poly([1]),m-2*X]
        for k in range(1,len(coeff)-1):
            ls.append(((2*k+m-2*X)*ls[-1]-(k+m-1)*ls[-2])/Q(k+1))
        self.P = sum((Q(str(F(str(c))))*ls[k]/ls[k][0] for k,c in enumerate(coeff)),Poly([]))
        self.P1 = self.P.derivative()-self.P
        self.P2 = self.P1.derivative()-self.P1
        self.pi = arb.pi()

    def f(self,d):
        x = self.pi*arb(d)**2
        return ev(self.P,x)*(-x).exp()

    def numerical_critical_radii(self):
        result = []
        for z,multiplicity in self.P1.complex_roots():
            if z.imag.is_zero() and z.real > 0:
                result.append(float((z.real/self.pi).sqrt()))
        return sorted(result)

    def verify(self,lo,hi,B,eta,max_depth=45,max_cells=200000):
        """Prove f(c)+B*eta<0 throughout the CLOSED rational interval."""
        stack = [(lo,hi,0)]
        checked = maximum_depth = 0
        failure = None
        while stack:
            a,b,depth = stack.pop()
            center = (a+b)/2; h = (b-a)/2
            ca = arb(Q(str(center))); ha = arb(Q(str(h)))
            x = self.pi*ca*ca; ex = (-x).exp()
            g = ev(self.P,x)*ex+arb(Q(B)*Q(str(eta)))
            gp = 2*self.pi*ca*ev(self.P1,x)*ex
            xm = self.pi*arb(Q(str((a*a+b*b)/2)))
            xh = self.pi*arb(Q(str((b*b-a*a)/2)))
            M2 = (2*self.pi*mag(self.P1,xm,xh)
                  + 4*self.pi*self.pi*arb(Q(str(b*b)))*mag(self.P2,xm,xh))
            M2 *= (-self.pi*arb(Q(str(a*a)))).exp()
            upper = g+abs(gp)*ha+M2*ha*ha/2
            checked += 1; maximum_depth = max(maximum_depth,depth)
            if upper < 0:
                continue
            if g > 0 or depth >= max_depth or checked >= max_cells:
                failure = dict(lo=str(a),hi=str(b),center_value=str(g),upper=str(upper),depth=depth)
                break
            stack.extend([(a,center,depth+1),(center,b,depth+1)])
        return dict(verified=failure is None and not stack,checked_cells=checked,
                    max_depth=maximum_depth,failure=failure)


def grid(width,cap_geometry=None):
    if cap_geometry is None:cap_geometry=geometry()
    lo = F(cap_geometry[0]['lower_base']); hi = 2*RADIUS
    assert width > 0
    points = {lo,hi,*[F(row['lower_base'])for row in cap_geometry]}
    p = lo
    while p < hi:
        points.add(p)
        p += width
    return sorted(points)


def construct(n,width,safety,mode='centroid'):
    start = time.monotonic()
    f = PhysicalFunction(SOURCE/SOURCE_NAMES[n])
    critical = f.numerical_critical_radii()
    cap_geometry=geometry(n,mode)
    edges = grid(width,cap_geometry)
    bins = []
    for lo,hi in zip(edges[:-1],edges[1:]):
        B = cap_from_geometry(lo,cap_geometry)
        points = [float(lo),float(hi)] + [x for x in critical if float(lo)<x<float(hi)]
        # High-precision values avoid catastrophic binary64 cancellation near
        # the nearly double roots. Selection is heuristic; verification below
        # is the authority even if this critical-point selection were wrong.
        minimum = min(float(-f.f(arb(Q(str(F(str(x)))))))/B for x in points)
        assert minimum > 0,(n,lo,hi,minimum)
        eta = F(format(float(safety)*minimum,'.12g'))
        attempts = []
        for attempt in range(12):
            check = f.verify(lo,hi,B,eta)
            attempts.append(dict(eta=str(eta),**check))
            if check['verified']:
                break
            eta /= 2
        assert check['verified'],(n,lo,hi,attempts)
        bins.append(dict(lo=str(lo),hi=str(hi),cap=B,eta=str(eta),
                         numerical_minimum_budget=minimum,
                         approximate_budget_fraction=float(eta)/minimum,
                         verification=check,construction_attempts=len(attempts)))
    boundaries = []
    for i,d in enumerate(edges):
        left = F(bins[i-1]['eta']) if i else F(0)
        right = F(bins[i]['eta']) if i<len(bins) else F(0)
        boundaries.append(dict(base=str(d),eta=str(min(left,right))))
    return dict(dimension=n,source=str(f.path.relative_to(ROOT)),source_sha256=f.sha256,
                physical_function=f.payload,physical_function_sha256=f.physical_sha256,
                short_radius=str(RADIUS),width=str(width),construction_safety_fraction=str(safety),
                geometry_mode=mode,
                precision_bits=ctx.prec,verified=True,bins=bins,boundaries=boundaries,
                exact_cap_geometry=cap_geometry,seconds=time.monotonic()-start,
                convention='Use bin eta only on its open interior. At each endpoint use boundaries eta; outside support use zero.',
                verifier='Exact rational Laguerre reconstruction; 256-bit Arb second-order Taylor bounds on closed bins.',
                scope='Pair-budget and geometric-cap certificate only. No new kernel or density bound is certified.')


def replay(path,candidate=None):
    data = json.loads(path.read_text())
    records = data['dimensions']
    report = []
    for item in records:
        # Reconstruct directly from the portable rational payload. Full-source
        # SHA is retained as provenance, not as the budget's applicability key.
        f = PhysicalFunction(payload=item['physical_function'])
        assert f.n==item['dimension'] and f.physical_sha256==item['physical_function_sha256']
        assert F(item['short_radius']) == RADIUS
        mode=item.get('geometry_mode','centroid')
        codes=[row['spherical_certificate']for row in item['exact_cap_geometry']if 'spherical_certificate'in row]
        cap_geometry=geometry(f.n,mode,codes)
        assert item['exact_cap_geometry'] == cap_geometry
        edges = [F(item['bins'][0]['lo'])]
        checks = []
        for row in item['bins']:
            lo,hi,eta = F(row['lo']),F(row['hi']),F(row['eta'])
            assert lo == edges[-1] and lo<hi and eta>=0 and row['cap']==cap_from_geometry(lo,cap_geometry)
            assert not any(lo<F(row['lower_base'])<hi for row in cap_geometry)
            check = f.verify(lo,hi,row['cap'],eta)
            assert check['verified'],(item['dimension'],row,check)
            checks.append(check); edges.append(hi)
        assert edges[0]==F(cap_geometry[0]['lower_base']) and edges[-1]==2*RADIUS
        assert len(item['boundaries'])==len(edges)
        for i,(edge,boundary) in enumerate(zip(edges,item['boundaries'])):
            left=F(item['bins'][i-1]['eta']) if i else F(0)
            right=F(item['bins'][i]['eta']) if i<len(item['bins']) else F(0)
            assert F(boundary['base'])==edge and F(boundary['eta'])==min(left,right)
        report.append(dict(dimension=item['dimension'],verified=True,bins=len(checks),
                           checked_cells=sum(c['checked_cells']for c in checks),
                           max_depth=max(c['max_depth']for c in checks)))
    out=dict(verified=True,source=str(path),results=report)
    if candidate is not None:
        cf=PhysicalFunction(candidate)
        row=next(r for r in records if r['dimension']==cf.n)
        assert cf.physical_sha256==row['physical_function_sha256'],'Candidate physical function does not match budget'
        out['candidate_binding']=dict(candidate=str(candidate),dimension=cf.n,
                                      physical_function_sha256=cf.physical_sha256,verified=True)
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--width',default='1/200')
    ap.add_argument('--safety',default='17/20')
    ap.add_argument('--dimensions',type=int,nargs='+',default=[5,6,7])
    ap.add_argument('--output',type=Path,default=HERE/'step_allowances_005.json')
    ap.add_argument('--verify',type=Path)
    ap.add_argument('--candidate',type=Path,help='Optional candidate to bind by exact physical coefficients, independent of kernel matrices')
    ap.add_argument('--geometry',choices=['centroid','midpoint_rankin','transverse_rankin','spherical_transverse'],default='centroid')
    args=ap.parse_args()
    if args.verify:
        out=replay(args.verify,args.candidate)
        print(json.dumps(out,indent=2))
        args.verify.with_name(args.verify.stem+'_replay.json').write_text(json.dumps(out,indent=2)+'\n')
        return
    width=F(args.width);safety=F(args.safety)
    assert 0<safety<1
    rows=[]
    for n in args.dimensions:
        row=construct(n,width,safety,args.geometry);rows.append(row)
        print(dict(dimension=n,bins=len(row['bins']),seconds=row['seconds'],
                   min_eta=min(float(F(b['eta']))for b in row['bins']),
                   max_eta=max(float(F(b['eta']))for b in row['bins']),
                   max_depth=max(b['verification']['max_depth']for b in row['bins']),
                   fallback_bins=sum(b['construction_attempts']>1 for b in row['bins'])),flush=True)
    result=dict(schema='tail_funded_lens_step_allowances_v1',dimensions=rows)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Wrote',args.output)


if __name__=='__main__':main()
