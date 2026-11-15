# Resolved comparisons and further lower bounds, 15 November 2026

**Current literature status:** our 84D, 85D and 96D bounds are best-known
in the review dated 15 November 2026. No stronger published bound was located.
[Sources, exact comparisons and criterion](../covering_and_deep_cosets_20261115/PRIORITY.md).
Older checkpoint and candidate labels below describe the earlier investigation.

**Current update:** the [covering-based extension](../covering_and_deep_cosets_20261115/)
raises the 96D density to **1.039108834467248e16**, or **118.554980%**
above BCS1982. The original 96D figures below are a superseded checkpoint;
their proof and source reconstruction remain available.

The subsequent [complete 48D–96D audit](../range48_96_20261115/)
adds explicit 49D–52D lattices, a stronger 81D bound and a sourced table
for every dimension. The 84D/85D/96D comparisons below are retained.

**Original 96D checkpoint:** a periodic packing exists with center density
**8,236,644,396,376,130.309…**, **73.240722% above** the published
Bos–Conway–Sloane 96D construction. This is an existence proof with nine
cosets; explicit horizontal offsets have not been found.

The same investigation corrects the earlier comparisons. The classical
BCS88 construction has sections stronger than our 86D and 87D results,
so those are **not record candidates**. Our 84D and 85D bounds still beat
these stronger classical comparisons, by **156%** and **32.197832%**.
The earlier larger percentages against Chen remain arithmetically valid
but use a weaker benchmark.

All numbers below are **center densities**. Multiply by
`V_n = pi^(n/2)/Gamma(1+n/2)` to obtain occupied fractions.

| Dimension | External comparison | Our bound | Verdict |
| --- | ---: | ---: | --- |
| 84 | `2^39 = 549,755,813,888`, BCS88 section, **derived** | `2^45/25 = 1,407,374,883,553.28` | **+156%**; record candidate against the strongest comparison located |
| 85 | `2^40 = 1,099,511,627,776`, BCS88 section, **derived** | `2^47/(25√15) = 1,453,530,529,565.988…` | **+32.197832%**; record candidate against the strongest comparison located |
| 86 | `2^42/√3 = 2,539,213,337,094.388…`, BCS88 section, **derived** | L80 counting: 1,504,719,014,574.453…; 88D antipode: 1,540,370,290,434.365… | **Both rejected as record candidates** |
| 87 | `2^(85/2) = 6,219,777,023,950.949…`, BCS88 section, **derived** | 88D antipode: 4,551,964,265,436.057… | **Rejected as a record candidate** |
| 88 | `2^44`, already obtained by BCS's classical theorem | Recent Sakamoto lattices also have `2^44` | No density improvement; full upstream minimum computation now reproduced |
| 93 | MW92 × A1: 4.887025440742e15, **derived product** | **8.182565090782e15** | New periodic existence bound; no reported dimension-specific record identified |
| 94 | MW92 × A2: 2.821525453749e15, **derived product** | **7.702212346015e15** | Same limited comparison scope |
| 95 | MW92 × A3: 1.727824414490e15, **derived product** | **7.410425802173e15** | Same limited comparison scope |
| 96 | **BCS1982: `3^48/2^24 = 4,754,450,504,593,402.735…`, reported construction** | **8,236,644,396,376,130.309…** | Historical **+73.240722%** checkpoint, superseded by the covering-based bound linked above |

The original “record candidate” labels are superseded by the dated
best-known claims linked above. In 84D/85D the comparator is derived
from a classical theorem; in 96D the paper explicitly reports the
construction. The weaker 86D/87D candidates remain rejected.
The research ceiling remains **96D**, including all parent inputs.

## Proofs and verification

* [Periodic 93–96D proof](TORUS_PROOF.md),
  [exact checker](verify_torus_extensions.py),
  [exact density squares and gains](torus93_96_results.json).
  This applies the greedy-torus/colored-layer method already used in the
  repository, with better small-dimensional colorings; it does not claim
  a new construction principle.
* [BCS88 section derivation](CLASSICAL_SECTIONS.md),
  [finite-field and comparison checker](verify_classical_sections.py),
  [full matrix reconstruction](construct_bcs88.py),
  [Gram data](bcs88_data.npz). These comparisons do **not** depend on
  Sakamoto's new 88D minimum theorem.
* [BCS96 reconstruction and local-optimality proof](BCS96_LOCAL_PROOF.md),
  [exact matrix auditor](audit_bcs96.py). This closes the direct local
  lattice-deformation route at BCS96; it says nothing about global
  optimality or the denser periodic construction.
* [Fresh 88D source replay](REPLAY88.md): all eight top-level stages
  passed in 1,697.53 seconds. The rebuilt L0 Gram equals our archived
  catalogue Gram exactly. This is a fresh execution of the author's
  computational proof, not a new implementation or a formal proof.

From the repository root, with assertions enabled:

```sh
python research/priority_resolution_20261115/verify_torus_extensions.py
python research/priority_resolution_20261115/verify_classical_sections.py
python research/priority_resolution_20261115/construct_bcs88.py
python research/priority_resolution_20261115/audit_bcs96.py
```

The first two use the standard library. The matrix auditors also need
NumPy and python-flint. Saved arithmetic reports are checked for exact
agreement; `--write-result` deliberately regenerates them.
Published parent-minimum theorems, the BCS construction theorem and the
classical root-lattice facts are stated mathematical inputs.

## Sources and search scope

1. H. J. A. Bos, J. H. Conway and N. J. A. Sloane,
   [*Further lattice packings in high dimensions*](https://neilsloane.com/doc/Me103.pdf),
   Mathematika **29** (1982), 171–180,
   [DOI](https://doi.org/10.1112/S0025579300012262).
   Theorems 1–2 and equation (15) give the construction; Table 1, p.177,
   reports the 96D `log2(delta)=52.078` entry from P48q; §5 gives its
   norm-doubler. Exact reconstruction gives `delta96=3^48/2^24`.
   The downloaded PDF SHA-256 is
   `8ce46a5a6af6e58bd2a2f7ab2ceefcd629a106b9c61b0ddae99c59940054d448`.
   The authors' determinant convention is covolume; our Gram determinant
   is its square.
2. M. Schütt and T. Shioda,
   [*Mordell–Weil Lattices*, author draft 16 March 2017](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf),
   Theorem 13.35, equation (13.24), p.357: at `p=47`, rank 92,
   minimum 16 and Gram determinant `47^14`. This is the MW92 parent,
   credited to Shioda, used in the new periodic proof.
3. R. Sakamoto, [arXiv:2609.15633v2](https://arxiv.org/abs/2609.15633v2),
   [source and reproducibility supplement](https://github.com/RyotaroSakamoto0113/rank-88-extremal-lattice-computations/tree/628ce2a7ac23acf8d4c295677e4310c9ec5b21d9).
   See the pinned replay report for the resolved minimum-computation gap.
4. The [earlier literature review](../cross_dimension_lower_20261115/LITERATURE_STATUS.md)
   records Chen, Shimada, code tables, Cohn, Nebe–Sloane and citation
   searches. This continuation additionally recovered the complete BCS
   paper from Sloane's publication archive and searched for later 96D
   improvements. The [four-Leech 96D page](https://srfp311t1.com/rank96_four_leech)
   explicitly lacks a complete minimum-10 certificate, so its candidate
   is not adopted as a theorem. No stronger proved numerical 96D bound
   was located. This search result supports the candidate label, rather
   than proving absence of every competing construction.

## Next mathematical targets

Find explicit torus offsets for the 96D nine-color proof, or exploit
MW92's structure to obtain deeper compatible offsets. In 84D/85D,
produce explicit representatives for the existing binary-counting
proofs and compare against further structured MW92 sections. In 86D/87D,
any renewed search must first beat the classical numbers above. Direct
infinitesimal lattice deformation of BCS96 cannot produce a gain.
