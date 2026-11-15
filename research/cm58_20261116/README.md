# A new explicit best-known lower bound in 58 dimensions

**Center density: 161,833.931161289966…**, an improvement of
**6.927715744…%** over the stronger classical Gamma72/Leech-block
section, whose center density is 151,348.908966445256… .

The comparison is derived from Nebe's published 72D lattice and a
primitive rank 14 Leech-block normal. Its comparison formula and original source are recorded explicitly.
The new lattice also improves the directly reported Elkies 58D guarantee
in [Vardy (1999), p.278](https://doi.org/10.1090/S0002-9947-99-02169-8)
by **>78.281202%**, using the upper end of its published logarithmic prefix.
No stronger published 58D guarantee was located in the dated review.

The construction is an **explicit lattice with minimum exactly 10**.
Its [complete rational 58D Gram](gram58.json) and [certificate](certificate.npz)
are supplied. The parent is Shioda's published 56D Mordell-Weil lattice;
its minimum 10 is inherited, while the explicit source-to-Gram binding
is checked from 120 polynomial sections and 7,140 pairwise heights.

The useful shortcut is an order-three isometry C with C²+C+I=0.
One certified deep half-coset supplies all three nonzero classes of a
compatible binary plane. A scaled A2 vertical lattice completes the
58D construction. Two complete enumeration runs each visited
**2,904,500,151 nodes**. Both completed reports and the exact negative
control are archived. The replay uses the same audited implementation.

```sh
python research/cm58_20261116/verify.py
python research/cm58_20261116/verify.py --replay /tmp/new_cm58_replay
```

NumPy, SymPy and python-flint are required; the full replay also needs
GNU C++. Run without `-O`. The replay writes into a new directory and
runs two complete positive searches plus the negative control. The
saved-report audit is quick; the full searches take several minutes each.

[Proof and arithmetic assumptions](PROOF.md) · [Sources](sources.json) ·
[Exact result](results.json) · [The other nine research targets](FRONTIER.md).
