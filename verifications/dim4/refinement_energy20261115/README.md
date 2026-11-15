# 4D: a stricter bound than the Felpix certificate

**Every sphere packing in R4 has upper density < 0.629520189662.**
The pinned Felpix-Studios result is <0.629521681653. The reduction between
the strict bounds is **0.000001491991**. The D4 lower bound remains pi^2/16;
optimality is not established.

The new ingredient is a degree-16 spherical energy polynomial supported
by five unused angular Gram credits. Its full three-variable certificate
has **308 interval leaves**, and **245 exact Gram range equations** justify
the negative polynomial coefficients. The larger-radius contribution
retains the original certified coefficient. See [the proof](PROOF.md) and
[attribution](NOTICE.md).

## Reproduce

Use Python 3.12 or later, the pinned dependencies, and a C++20 compiler:

```sh
python -m pip install -r verifications/dim4/refinement_energy20261115/requirements.txt
CXX=g++ python verifications/dim4/refinement_energy20261115/verify.py \
  --out /tmp/d4-energy-full --workers 2
```

The output directory must be new. The default command performs a fresh
full upstream replay and the new checks. It requires no solver or network
after dependency installation. The complete run takes substantially longer
than the small extension. Do not use `python -O` or `PYTHONOPTIMIZE`.

The [release audit](audit/) records the fresh full base replay (about 38
minutes here), composition with its fresh census, exact export checks and
six rejected invalid inputs. Other dimensions were unchanged in this release.

To check just the extension, with its scope explicit:

```sh
python verifications/dim4/refinement_energy20261115/verify.py \
  --extension-only --out /tmp/d4-energy-extension
```

`--base-replay DIR` instead consumes a separately completed upstream full
replay. The output records reuse; it does not claim the base was freshly
executed by that invocation.

## Files

- `certificate.json`: exact final bound, binding gain and frozen input hashes.
- `polynomials.json`: exact normalized Gegenbauer coefficients, scale included.
- `gram_credits.json`: exact range solutions in the original Gram blocks.
- `kernel_box_j1.json`: complete new subdivision tree.
- `verify.py`, `gram.py`, `polynomial.py`, `boxes.py`: replay and production tools.
- `vendor/`: unchanged, attributed upstream proof and its license.

`boxes.py --j 1 --out NEW.json` regenerates the new interval tree. Numerical
search is separate from verification. Comparison is to the stated pinned
repository, with no exhaustive worldwide-priority assertion.
