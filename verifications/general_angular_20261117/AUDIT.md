# Independent review of unrestricted45D/48D angular substitutions

Status: PASS. Both compact angular gates, exact structure/tails, both full radial/Fourier paths and strict rounding checks completed successfully. 45D density <1.2613300683e−6; 48D density <4.0691468600e−7. No publication is performed by this package.

## Frozen scope

The angular SHA256 hashes are45D ec1212ebad3ffbfcaae6ccff54d7eb61158549c0b2f127cb5e935b60242288d9 and48D862b3b7a2c472ad666425e1bcff3e5a8d45810c0ac0a1fea54492642889882f2. Source P,H,s,a,b,q and radial cutoff are unchanged. Exact Cnew/Bnew=Cold/Bold, while C decreases. This leaves the entire pair inequality unchanged and decreases the origin contribution. check_substitution.py proves these identities and final density rounding using separately constructed201/200-term rational Machin bounds.

## Theorem and global coverage

Directly reviewed Cohn–de Laat–Salmon, arXiv2206.15373, §2.2 equation(2.4), Theorem2.4. The cyclic three-point condition, not the stronger pointwise premise of Theorem1.1, is used. The source's harmonic feature construction yields a PSD kernel, including the origin block [[C,−C/B],[−C/B,C/B²]]. A smooth cutoff inside the excluded minimum-distance ball removes the angular singularity without changing any packing-domain value. Gaussian decay provides Schwartz regularity. Thus this is an unrestricted packing theorem, not a lattice-only argument.

For ordered triangle side lengths sqrt(sT)*(1,u,w), the full domain is T≥1,u≥1,−1≤t≤1/(2u), w²=u²+1−2ut. The compact T=1 domain is exactly the union of two closed rectangles: negative t in[−1,0],u in[1,8]; and positive t=v/(2u), v in[0,1],u in[1,8]. The transformed gap1−v is algebraically exact. Boundaries and degenerate collinear triangles are included, not discarded.

At T=1 the verified expression is Q(t)+A+B. At general T the secondary terms acquire factors exp[−216(T−1)g1] and exp[−216(T−1)g2], with0≤g1≤g2. The secondary angles satisfy beta≥alpha because (u−1)(1+t)≥0. The exact unique negative-to-positive sign transition excludes A>0,B<0. When A≤0≤B, the damped sum is at most max(0,A+B); when both positive it decreases, and when both negative Q(t)<0 suffices. Thus compact T=1 verification plus the strict negative cap proves all dilations, not just sampled scales.

For u≥8, Q(t)≤−delta. A positive alpha term requires alpha>cut, hence gap>2*cut*u−1; w/u≤2 gives the first tail bound4^192 exp[−216(2*cut*u−1)]. The other term is bounded by u^384 exp[−216(u²−1)], decreasing on this range. Multiplication by B bounds angular magnitudes using nonnegative normalized Gegenbauer coefficients. The exact structure verifier proves this combined tail<delta/2 and verifies U=8. This covers every shape outside the compact domain.

## Interval implementation review

Reviewed every Jet operation: product includes both mixed-gradient terms; reciprocal uses2g_i g_j/v³−h_ij/v²; logarithm uses h_ij/v−g_i g_j/v²; exponential uses exp(v)*(h_ij+g_i g_j). Polynomial composition uses the exact derivative chain rule, evaluated by Arb after shifting to the interval midpoint. All these give outward derivative enclosures.

For each rational rectangle, Taylor's theorem gives f(center)+sum|grad(center)_i|r_i+1/2 sum|H_ij(box)|r_i r_j. Summing all ordered Hessian indices includes the mixed term with its correct multiplicity. Direct interval and first-order enclosures are alternative sufficient tests. The split score is heuristic only; every midpoint is exact rational and either child partitions the closed parent. The audit copy rejects disabled assertions, verifies full polynomial degree, requires empty stack, and checks the two-root forest identity nodes=2*leaves−2. It has no success path at an iteration cap.

The reviewed producer compact implementation was replayed at768bits:45D9910nodes/4956leaves/maxdepth28;48D10146nodes/5074leaves/maxdepth32. Exact harmonic reconstruction, cap, Sturm sign transition and tail checks were replayed at1024bits. This is an independent review and higher-precision replay of that implementation, not a separately authored compact solver or formal proof. Both arithmetic paths share FLINT/Arb.

The radial checks use separate direct Fourier-coefficient and Laguerre-recurrence reconstructions, complete polynomial root accounting for global rays, full interval pair coverage, origin PSD and exact decimal rounding. Audit copies use1024bits. The original legacy angular million-box programs were not completed and are not claimed to have passed; the reviewed standalone compact proof and exact structure/tail gate replace them.
