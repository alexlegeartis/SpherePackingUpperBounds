# CT36 contact-preserving symmetry-breaking flex

2026-11-12. **UNVERIFIED packing candidates, not density records.**

## New mechanism

Instead of changing only three scalar eigenspace scales, allow the full symmetric Gram matrix to move while preserving its known contacts. This is related in spirit to the contact-accumulating ellipsoid deformation in Klartag, *Lattice packing of spheres in high dimensions using a stochastically evolving ellipsoid*, Inventiones244(2026),1251–1279, https://doi.org/10.1007/s00222-026-01412-w . Its asymptotic theorem does not itself imply a finite-dimensional record here.

Starting from the CTturyn2 three-scale vertex(a,b,c)=(1/2,1,5/2), the saved35864 norm8 vectors constrain the666 Gram coefficients with numerical rank665. A numerical null vector reconstructed to a small integer symmetric matrix D. **Every saved contact satisfies v^T D v=0 exactly.** This direction is outside the previously excluded three-scalar family.

Base source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CTturyn2.html . Base transformation and incomplete shell provenance are in sibling ct36_metric_screen_20261112 and ct36_symmetry_break_20261112. The35864-vector list came from a capped floating enumeration, not a complete shell certificate.

## Exact matrix facts, independently audited

Let G=G4/4 from flex.npz and G(t)=G+tD. The sibling audit_ct36_flex independently verified:

- detG=3^18(5/4)^12 and positive-definite exact LDL;
- A=G^-1D satisfies A³=(16/15)A, rank24, trace0;
- characteristic polynomial λ^12(λ²−16/15)^12;
- therefore detG(t)=detG(1−16t²/15)^12;
- G(t) is positive definite for |t|<sqrt15/4.

If its GLOBAL minimum remains8, its density/current36D-record ratio is

    (12288/15625) / (1−16t²/15)^6.

The improvement threshold is|t|>0.19182618074. At t=1/5 the ratio is exactly2187000000000000/2140753641621841≈1.02160283999: conditional **+2.160284%**. Candidate integer Gram20G(1/5)=5G4+4D is candidate_1.npy. candidate_-1.npy is the opposite direction.

At t=1/4 the conditional ratio is approximately1.18970837, or+18.97%. The independent audit tested90million explicit neighbors V±e_i±e_j and found a contact event at|t|=1/4: an oldnorm17/2 vector with slope−2 reaches8. It does not establish that no other vector intervenes sooner. See audit_ct36_flex/pairs.json.

## Bounded searches, not minimum proofs

- Both t=±1/5 tested in the original and a furtherLLL-reduced basis, with200million-node caps per run; no vector below8 found before any cap.
- Sixteen independently randomized unimodular/LLL bases at t=1/5,15million-node cap each; again no rejection.
- All searches use floating Cholesky pruning and are incomplete. Neither absence of a discovered short vector nor preservation of the partial contact list certifies minimum8.
- Matrix and contact facts are exact; the global minimum is the sole outstanding packing-certification gate.

## Files

flex.npz: exact G4,D,V. flex.py: discovery plus exact contact residuals. reduced_±1.npz: exact basis changes and integer candidate Grams. random_bases.json: all16 unimodular transforms and capped-search outcomes. short.cpp: bounded witness search; a discovered integer witness is evaluated exactly, but a capped run is not proof. The independent audit owns exact determinant and contact proofs in research/audit_ct36_flex/.
