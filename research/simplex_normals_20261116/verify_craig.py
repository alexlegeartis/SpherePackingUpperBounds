"""Exact ceilings on Chen's binary-code density GUARANTEE formula.

These are not upper bounds on actual Craig-lattice packing densities.
LP certificates are checked as rational Krawtchouk polynomials; no LP solver
or floating arithmetic is used. The 92D code case is explicitly unresolved.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
import argparse,json

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
DIMENSIONS = [48,54,55,56,63,64,72,80,92]


def prime(n):
    return n>=2 and all(n%j for j in range(2,isqrt(n)+1))


def kw(n,j,w):
    return sum((-1)**t*comb(w,t)*comb(n-w,j-t)
               for t in range(max(0,j-(n-w)),min(j,w)+1))


def polynomial_bound(cert,N,d):
    shortened = cert['shortened']; n = cert['length_after_shortening']
    assert isinstance(shortened,int) and 0<=shortened<N and n == N-shortened
    degrees = cert['krawtchouk_degrees']; coefficients = list(map(F,cert['normalized_coefficients']))
    assert len(degrees) == len(coefficients) and degrees[0] == 0 and coefficients[0] == 1
    assert degrees == sorted(set(degrees)) and all(0<=j<=n for j in degrees)
    assert all(a>=0 for a in coefficients)
    # An even binary code has only even pairwise distances. Shortening at
    # zero preserves evenness, minimum distance, and at least k-r dimensions.
    for w in range(d,n+1,2):
        assert sum(a*F(kw(n,j,w),comb(n,j))for j,a in zip(degrees,coefficients))<=0
    bound = 2**shortened*sum(coefficients)
    assert bound == F(cert['exact_code_size_bound'])
    return bound


def density_square(n,m,ell,k):
    return F(2)**(2*k-n)*m**n/F(ell**(2*m-2)*(n+1))


def audit():
    benchmarks = {int(n):F(x)for n,x in json.loads((P/'craig_benchmarks.json').read_text())['dimensions'].items()}
    assert sorted(benchmarks) == DIMENSIONS
    expected = {48:F(3,2)**48,54:F(5,2)**54/3**25,
                55:F(5,2)**54/3**25*F(9,4),56:F(5,2)**56/29**8,
                64:F(3**32),72:F(2**72),80:F(7,2)**80/41**12,
                92:F(2**184,47**14)}
    assert all(benchmarks[n] == v for n,v in expected.items())
    floor63 = benchmarks[63]
    assert floor63.denominator == 1 and floor63**25 <= 2**1196 < (floor63+1)**25
    data = json.loads((P/'code_certificates.json').read_text())
    assert [(a['dimension'],a['m'])for a in data] == [(n,m)for n in DIMENSIONS for m in range(1,(n+1)//2+1)]
    results = []; unresolved = []; polynomial_checks = 0
    for n in DIMENSIONS:
        maximum = F(0); maximizer = None
        for a in [x for x in data if x['dimension'] == n]:
            N = n+1; m = a['m']; d = 8*m
            ell = next(x for x in range(N,2*N)if prime(x))
            assert a['binary_length'] == N and a['distance'] == d and a['ell'] == ell
            # Larger admissible primes only reduce the formula (or leave
            # it unchanged at m=1), so the smallest prime covers them all.
            g = max([0]+[k for k in range(1,n+1)
                         if sum((d+(1<<j)-1)//(1<<j)for j in range(k))<=N])
            volume = sum(comb(N,j)for j in range(min((d-1)//2,N)+1))
            h = max(0,N-(volume-1).bit_length())
            assert g == a['griesmer_dimension_upper'] and h == a['hamming_dimension_upper']
            kmax = min(g,h) if d<=N else 0
            assert kmax == a['simple_dimension_upper']
            if a.get('LP_certificate') is not None:
                bound = polynomial_bound(a['LP_certificate'],N,d); polynomial_checks+=1
                lp = max(k for k in range(N+1)if 2**k<=bound)
                kmax = min(kmax,lp)
            candidate = density_square(n,m,ell,kmax)
            if candidate > maximum:
                maximum,maximizer = candidate,{'m':m,'ell':ell,'code_dimension_upper':kmax}
            if candidate > benchmarks[n]:
                needed = next(k for k in range(n+1)if density_square(n,m,ell,k)>benchmarks[n])
                unresolved.append({'dimension':n,'m':m,'ell':ell,'even_binary_length':N,
                                   'minimum_distance_required':d,'dimension_needed':needed,
                                   'dimension_upper':kmax,'conditional_density_squared':str(candidate),
                                   'conditional_ratio_squared':str(candidate/benchmarks[n])})
        results.append({'dimension':n,'guarantee_formula_density_squared_ceiling':str(maximum),
                        'comparison_density_squared':str(benchmarks[n]),'maximizing_parameters':maximizer,
                        'formula_cannot_improve_comparison':maximum<=benchmarks[n]})
    assert len(unresolved) == 1
    u = unresolved[0]
    assert (u['dimension'],u['m'],u['ell'],u['even_binary_length'],u['dimension_needed'],u['dimension_upper'],u['minimum_distance_required']) == (92,3,97,93,43,43,24)
    assert [a['dimension']for a in results if a['formula_cannot_improve_comparison']] == DIMENSIONS[:-1]
    return {'status':'eight_formula_families_ruled_out_one_specific_code_case_open',
            'scope':'Ceilings on the guarantee formula in Chen2013 Theorems2.3-2.4, using even binary linear codes; not ceilings on actual lattice densities or other constructions.',
            'rows':results,'polynomial_certificates_checked':polynomial_checks,
            'all_polynomial_and_density_checks_exact':True,
            'unresolved_conditional_cases':unresolved,'new_packing_SOTA_claimed':False,
            'maximum_parent_or_target_dimension':93}


if __name__ == '__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--write-result',action='store_true'); args=ap.parse_args()
    result=audit()
    if args.write_result: (P/'craig_results.json').write_text(json.dumps(result,indent=2)+'\n')
    else: assert json.loads((P/'craig_results.json').read_text()) == result
    print(result['status'],'rational LP certificates',result['polynomial_certificates_checked'])
