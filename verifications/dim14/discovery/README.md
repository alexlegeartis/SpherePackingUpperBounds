# 14D discovery provenance

The numerical program is pure degree-63 normalized Laguerre radial expansion,
with scalar angular degree at most 19. Minimum squared separation is 7/4.
Profiles (a,b)=(12,16),(14,8),(16,0) each received two angular separation rounds,
including the full negative-angle interval [-1,1/2]. The 15D reference under
inputs/ supplies positive numerical column scales only; all Gegenbauer
recurrences, constraints and optimization were recomputed in dimension 14.

The numerical screening densities were respectively 0.04500316036545258,
0.04499999689853886 and 0.04502411847821455. None is a bound. The first two
looked competitive against 0.04501031; the matched single-Gaussian control
was noncompetitive. The better curved profile (14,8) was frozen and selected
for exact proof. The other curved profile was not subjected to heavy proof
work; this is not an exhaustive profile optimization.

For the selected angular candidate, exact binary input coefficients receive
a constant guard 1/10000 and normalization to constant harmonic coefficient
one. The full angular proof then passes. The radial search uses that guarded
B, physical cutoff X=12, margin 1e-8 and the saved cut file. Every final scalar
coefficient is rationalized by Fraction(float(value)).limit_denominator(10^16).
The origin multiplier receives exactly the factor 1000001/1000000.

Initial sampled radial candidates failed exact Fourier or physical root tests.
Their roots and interval midpoints were added to the LP. Round 7 first passed
both rays, but a narrow positive pair lobe remained near X=9.474. After round
8, a validated point evaluation still found it positive. The actual local
maximum and nearby points were added in round 9; complete root and pair
proofs then passed. All rejections and the final success remain in
history_cut12/ and refinement_cut12_history.json.

Optional numerical replay needs numpy, scipy and python-flint. From this
folder, representative commands are:

- python screen.py --a 14 --b 8 --rounds 2
- python search_radial.py --dimension 14 --degree 63 --s 7/4 --a 14 --b 8 --cut 12 --margin 1e-8 --tag radial_cut12

Fresh solves may differ numerically. For acceptance of the frozen certificate,
use the parent folder's verify.py, which does not depend on optimizer output
being reproducible. Historical refinement scripts are provenance, not proof
entrypoints. The final candidate's angular path is local to this folder.
