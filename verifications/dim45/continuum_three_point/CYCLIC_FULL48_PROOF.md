# A complete cyclic continuum certificate in dimension 48

Producer verification gives normalized density objective 18.911756051853086... <18.912. Independent verification is requested. This is an instance of the existing Cohn–de Laat–Salmon continuum framework, not a new method or an optimality/record claim.

## Data

- Minimum center separation r=sqrt(6).
- Radial Fourier-positive function f2(x)=exp(-pi|x|²)P(pi|x|²), where P and its Fourier polynomial H are in adaptive_cyclic111.json, each divided by H(0). Degree111.
- Angular degree19 polynomial Q from angular_improved48_root501_exact.json, constant normalized Gegenbauer coefficient1, all coefficients nonnegative, B=Q(1).
- k=6 and h(x)=exp(-6pi(|x|²-6)) for |x|>=r.
- M is the exact rational neighbor_cost/(B H(0)), multiplied by1000001/1000000; its full value is in cyclic_radial111_exact.json.

## Positive definite kernel

Take a smooth radial cutoff chi, zero near the origin and one for |x|>=r. Put h(x)=exp(-6pi(|x|²-6)) globally. Set g0(x)=chi(x)h(x)-(1-chi(x))B and g(x)=chi(x)h(x). For nonzero x,y, define

 f3(x,y)=(M/B)[g0(x)g0(y)+g(x)g(y) sum_{l>=1} b_l P_l^(48)(xhat dot yhat)],

with the angular features zero at the origin. The addition theorem proves positive definiteness. All features are smooth and rapidly decreasing; the angular terms vanish in a neighborhood of the origin. Outside radius r, f3=(M/B)h(x)h(y)Q(angle).

At the required degenerate configurations,

 f3(0,0)=MB,
 f3(x,0)=f3(0,x)=-M h(x),
 f3(x,x)=M h(x)^2 for |x|>=r.

## Global cyclic sign

Rescale all triangle sides by r. The Gaussian parameter becomes alpha=36pi>96. CYCLIC96_PROOF.md and cyclic_triangle96_exact.json establish the required global cyclic sum is nonpositive for alpha96 and every larger alpha. Multiplication by the positive factor (M/B) exp(2alpha) does not change the sign.

## Radial conditions

The Fourier polynomial H is nonnegative on the whole nonnegative ray and H(0)=1, so fhat2>=0 with fhat2(0)=1. The certificate verifies that P is negative for pi|x|²>=7pi. This uses complete certified complex root isolation and exact endpoint signs, not sampling.

For 6<=|x|²<=7 the required coupled inequality is

 f2(x) <= M(2h(x)-h(x)^2).

Writing s=pi|x|², divide by h(x)>0. It suffices to check

 P(s) exp(5s-36pi) <= M(2-exp(-6(s-6pi))).

certify_cyclic_radial111.py verifies this inequality on the slightly larger rational interval [6*pi_lower,7*pi_upper]. Exact rational Machin bounds enclose pi. Each subdivision evaluates a centered polynomial Taylor expansion with Arb interval coefficients, bounding P from above by the constant coefficient plus the sum of absolute higher coefficients. Exponential monotonicity provides upper and lower bounds for the remaining factors. All1,931 intervals pass, with maximum depth21 and384-bit precision. For |x|²>=7, f2<=0 while 0<h<=1, so the same coupled inequality follows automatically.

## Consequence and scope

The strengthened continuum three-point theorem, Cohn–de Laat–Salmon Theorem2.4 and equation2.4 (https://arxiv.org/html/2206.15373), now bounds the packing density by vol(B_(sqrt(6)/2)^48)*(f2(0)+f3(0,0)). The exact normalized objective P(0)+BM is18.911756051853086... <18.912. Relative to the previous independently audited19.4739676269 certificate this is about2.89% smaller, subject to independent audit of this new certificate.

The method and theorem are existing literature. This is not a proof of optimal packing, an explicit packing construction, a kissing-number bound, a statement about192D, or a claim to have beaten the unknown optimal CE value. The public table comparison is separate and does not establish priority.
