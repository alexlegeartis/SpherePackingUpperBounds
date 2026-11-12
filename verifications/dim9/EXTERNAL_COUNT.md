# Explicit external theorem dependency

The 9D radial LP uses average ordered-neighbor count373 in the half-open shell
[1,1.003012), after normalizing minimum center separation to1. This shell is
strictly contained in the published [1,1.003013) shell.

Source: D. de Laat, F. M. de Oliveira Filho and F. Vallentin, Upper bounds for
packings of spheres of several radii, arXiv:1206.2608, §6/Table2.
https://arxiv.org/html/1206.2608#S6.T2

The paper states its bounds were rigorously checked. The original SDP is not
replayed by this reconstruction. This is an average degree bound, not a
pointwise kissing bound, and not unordered edge count per center. There is
no integer-floor operation. The objective term is373*lambda, with lambda>=0.
At the outer threshold, the no-count pair inequality is required.

The computational templates are not alone a full theorem proof: exact
angular PSD, all-scale cyclic reduction, full compact geometry coverage,
radial Fourier signs and pair inequalities, rounding and theorem linkage
must all be checked on the same frozen inputs before publication.
