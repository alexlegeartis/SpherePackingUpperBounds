# An explicit 53D section from Leterrier's published Gram

Leterrier's 2024 thesis gives an explicit rank 54 Gram G in Computational
Proposition 3.3.10, pp.168–169. The parent has minimum 10 and determinant
`3^25`. Theorem 3.4.1 supplies the triviality condition used for the stated
identification. These published parent and computational-Gram inputs are
inherited; we do not re-enumerate the parent's minimum vectors.

The displayed matrix was transcribed and checked against the downloaded
PDF. [Source locators, PDF hash and transcription hash](sources.json) ·
[Archived numerical matrix](sources/leterrier54_gram.txt).

The certificate supplies a primitive integer label vector a with

```
gcd(a_i)=1,
a*G^-1*a^T=10/3.
```

Its perpendicular period consists exactly of integer coefficient rows z
with `z*a^T=0`. An exact unimodular HNF transformation constructs a full
basis N of that kernel, rather than a finite-index sublattice. The
period Gram has

```
det(N*G*N^T)=det(G)*(10/3)=10*3^24.
```

It inherits minimum at least 10. The certificate also gives a vector in
the kernel of norm exactly 10, so this minimum is attained. Therefore

```
delta53 = (10/4)^(53/2)/sqrt(10*3^24)
        = 5^26/(2^27*3^12)
        = 20,890.80489885343… .
```

This is an explicit lattice. It replaces the earlier conservative
Hermite-only section guarantee `17,624.6617556196…`, improving that value
by about 18.53%. It is a new numerical section and certificate here; the
section principle and the 54D parent are classical/published inputs.

[Certificate](leterrier53.npz) · [Exact audit](verify_section53.py) ·
[Full period and results](section53_results.json).
The optional `--source-pdf FILE` also verifies the matrix extraction from
the pinned thesis PDF. A short-vector discovery run was capped; no claim
that the selected dual normal is a shortest dual vector is required.
