# Unrestricted compact cyclic harmonic-matrix family

This is a dimension-parametric admissible family, not a proved new numerical bound. Let n>=3, squared minimum distance1, R=19/13, h(t)=((R-t)/(R-1))_+^12 and H(t)=(h(t),h(t)^2). For every harmonic degree k use a PSD2x2 matrix B_k and positive conditioning coefficient c_k, with c_0=1 and q=sum c_k. Set

K(x,y)=q^-2 sum_k c_k P_k^(n)(cos(x,y)) H(|x|²)^T B_k H(|y|²),
K(0,x)=q^-1 z^T H(|x|²), K(0,0)=s.

The block [[s,z^T],[z,B_0]] PSD together with B_k PSD is a sufficient kernel PSD condition. Define Q(t)=sum c_k P_k^(n)(t) B_k. No parity/even-harmonic restriction is imposed.

For a genuine triangle with sorted squared side lengths1<=a<=b<=c, if b>=R all three cyclic terms vanish. If b<R<=c only the term using sides a,b survives. Its opposite cosine is (a+b-c)/(2sqrt(ab))<=(a+b-R)/(2sqrt(ab))<=1/2, using a,b<=R. It can be negative, so its sign must hold over the ENTIRE interval[-1,1/2], not only[0,1/2]. Sufficient conditions are
Q00(t)<=0,
Q00(t)+Q01(t)<=0,
Q00(t)+2Q01(t)+Q11(t)<=0.
These are the three distinct corner values of the bilinear expression Q00+(u+v)Q01+uvQ11 for u,v in[0,1], where u=h(a),v=h(b). Every actual pair lies in that rectangle. Boundary h=0 causes no division issue in this direct statement.

Only1<=a<=b<=c<=R remains. This compact cube consists of valid triangles because sqrt(R)<2. The cyclic condition is imposed directly using all three vertex angles and H(a),H(b),H(c). In the discovery model it is divided by h(a)h(b)>0 on interior samples, with separately constrained cutoff cases. A final proof must use a globally continuous formulation or rigorously handle the division/zero boundaries.

Let f(v)=e^-x P(x), x=pi|v|², with Fourier polynomial nonnegative everywhere and Fourier(f)(0)=1, contact x=d. In these coordinates use H(x/d). Writing W=sum c_k B_k/q, the radial residual is
f(v)+(H^T W H+2z^T H)/q.
There is NO lattice factor2 on the quadratic term and NO mod2 or antipodal budget. A full nonpositive residual gives number density<=P(0)+s. Certified nonnegative annular maxima g_i may instead be paid with the previously established universal per-center code caps, adding E=sum C_i g_i. Packing density then is (d/4)^(n/2)/Gamma(n/2+1)*(P(0)+s+E).

Only f needs the usual Schwartz/Poisson assumptions. The compact kernel enters finite local sums in a periodic packing; PSD and cyclic sums give the counting bound directly. Hard support is therefore not being used to truncate f or its Fourier sign domain.

Current screen.py enforces ONLY sampled compact geometry/cap/radial/Fourier conditions. It is a discovery optimizer. Promotion requires exact PSD, whole-interval corner polynomial proofs, complete compact cyclic-domain proof, exact Fourier positivity, all radial/tail signs and density normalization, plus independent audit. No certificate is claimed by this note.
