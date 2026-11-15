"""Render provenance-aware comparisons without presenting our bounds as literature."""
from pathlib import Path
import math
import sympy as sp


def link(label, target):
    return f"[{label}]({target})"


def render(report):
    sources = report["sources"]
    table = {row["dimension"]: row for row in report["rows"]}
    lines = [
        "# Published references and repository comparisons, dimensions 57–71",
        "",
        "**The strongest 59D comparison calculated here is now the classical kappa13 section "
        "of Gamma72, with center density 269064.727051…; our increase is 2.028457%. "
        "This comparison is derived from published inputs, not a 59D result quoted from Nebe's paper. "
        "The current literature review registers our 59D bound as best-known, with no stronger published bound located.**",
        "",
        "[Dated review and current comparisons](../../covering_and_deep_cosets_20261115/PRIORITY.md). "
        "Comparison is between proved lower bounds; a previous-record table and an older lattice's "
        "unknown actual minimum are not prerequisites.",
        "",
        "The same distinction applies to the earlier 65–71D percentage gains: their comparators "
        "were our derived constructions or an inherited-minimum section bound. "
        "They were not identified published same-dimensional records.",
        "",
        "All densities are **center densities** `delta`, not occupied-volume fractions. "
        "The latter are `Delta_n = pi^(n/2)/Gamma(n/2+1) * delta`. "
        "A comparison of exact rational squared densities proves the percentages below.",
        "",
        "## Which previous work is being used?",
        "",
        "| Reference | Actual role in this report |",
        "|---|---|",
    ]
    for sid in ("V99", "N98", "N99", "N12", "K13", "S19", "CS96"):
        source = sources[sid]
        title = f"{source['authors']}, *{source['title']}*"
        lines.append(f"| {link(sid + ': ' + title, source['url'])} | {source['role']} |")
    lines += [
        "",
        "Full publication details and precise locators are listed under **Source details** below. "
        "Parent-lattice papers and construction methods are credited as inputs; using a new "
        "section of a published parent does not mean we have improved that paper's density "
        "in its original dimension.",
        "",
        "## Published numerical comparison: Elkies's sections, reported by Vardy",
        "",
        f"{link('Vardy (1999), Section 4, p. 278', sources['V99']['url'])} reports lower bounds for "
        "Elkies's sections `E_n`, in a table headed `log2(delta(E_n))`. "
        "Here both numerical columns use **log2(center density)**, so their units agree. "
        "The ellipses are retained: the printed entries are not exact density formulas.",
        "",
        "| Dimension | Published `log2(delta(E_n))`, Vardy §4 | Strongest bound calculated here, `log2(delta)` | What the comparison establishes |",
        "|---|---:|---:|---|",
    ]
    for n in range(57, 64):
        row = table[n]
        ref, result = row["computed_reference"], row["repository_construction"]
        best = max(sp.Rational(ref["center_density_squared"]), sp.Rational(result["center_density_squared"]))
        label = row["published_reference"]["printed_log2_center_density"]
        comparison = "Below even the printed published lower bound; no record gain."
        if n != 63:
            comparison = "Improves the printed published lower-bound guarantee."
        if n == 59:
            comparison = "23.67–24.54% improvement over the published guarantee (verified interval)."
        lines.append(f"| {n} | {label} | {math.log2(float(sp.sqrt(best))):.6f} | {comparison} |")
    lines += [
        "",
        "There is **no percentage claim against the actual Elkies lattices** in this table. "
        "The original PDF was checked on 15 November 2026. Reference 10, p. 283 attributes "
        "the section results to **Elkies, personal communication, February 1996**. "
        "No exact E59 section parameters are supplied there.",
        "",
        "For 59D, let `B` be the guarantee reported by the printed prefix. Then "
        "`2^17.75 <= B < 2^17.76`, and exact integer-power comparisons prove "
        "**23.675981% < 100*(delta_ours/B-1) < 24.536216%**. "
        "The upper bound here is on the **reported guarantee**, not on the actual "
        "Elkies density. This establishes an improvement over the published lower bound; "
        "the unknown actual minimum is not required for that conclusion. "
        "[Source verification and derivation](../../../constructions/dim59/README.md#verified-external-comparison) · "
        "[Exact arithmetic report](../../../constructions/dim59/external_comparison.json).",
        "",
        "For **65–71D**, this audit has not identified a published same-dimensional numerical "
        "record to use as the comparator. Nebe's 64D and 72D parent results do not fill that gap. "
        "Their exact parent densities are respectively `3^16` and `2^36`; no improvement "
        "in either parent dimension is claimed here.",
        "",
        "## Exact comparisons with calculated reference constructions",
        "",
        "The 59D comparator is a consequence of **published Gamma72 structure and classical kappa13**. "
        "Its numerical 59D density is derived here. Other rows use the reference constructions "
        "calculated in this repository. The links identify each denominator; "
        "these percentages must not be described as gains over identified published records.",
        "",
        "| Dimension | Reference construction / source | Reference center density | Our certified construction | Center density | Increase |",
        "|---|---|---:|---|---:|---:|",
    ]
    for n in (59, 60, 61, 62, 65, 66, 67, 68, 69, 70, 71):
        row = table[n]
        ref, result = row["computed_reference"], row["repository_construction"]
        ratio2 = sp.Rational(result["center_density_squared"])/sp.Rational(ref["center_density_squared"])
        assert ratio2 > 1 and ref["origin"] == "computed_in_this_repository"
        percent = float(100*(sp.sqrt(ratio2)-1))
        ref_link = link(ref["name"], ref["evidence_path"])
        result_link = link(result["name"], result["certificate_path"])
        lines.append(f"| {n} | {ref_link} | {ref['center_density']:,.6f} | {result_link} | {result['center_density']:,.6f} | {percent:.6f}% |")
    lines += [
        "",
        "In particular, the 59D row is now the six-translate construction "
        "`274522.589741…`, replacing the obsolete bare-section value `244079.864699…` "
        "that the previous version called the current candidate. Its exact period, "
        "translate data and attained minimum are checked by "
        "[two implementations](../../../constructions/dim59/README.md#replay). "
        "They also verify the earlier internal reference `delta=2^18`. "
        "The stronger classical reference now used in this table is checked separately by "
        "[verify_literature_section.py](../../../constructions/dim59/verify_literature_section.py).",
        "",
        "The Gamma72 section references in 66–71D also attain a ceiling for the "
        "**bare-section density formula using inherited minimum 8**, via the known "
        "small-dimensional Hermite constants. That ceiling does not bound arbitrary "
        "packings or sections with a higher actual minimum, and it is not a literature record.",
        "",
        "## Controls that do not give record improvements",
        "",
        "| Dimension | Construction tested here | Center density | Why it is not an improvement |",
        "|---|---|---:|---|",
    ]
    reasons = {
        57: "Below our 56D-derived torus bound 151706.083320.",
        58: "Below our Gamma72 section 151348.908966.",
        63: "Below the Elkies lower bound printed in Vardy §4 (`log2(delta)=23.92…`).",
        64: "Below Nebe's published Ne64 density 43046721; this is a failed control, not our best 64D bound.",
    }
    for n in reasons:
        result = table[n]["repository_construction"]
        lines.append(f"| {n} | {link(result['name'], result['certificate_path'])} | {result['center_density']:,.6f} | {reasons[n]} |")
    lines += [
        "",
        "## Where our reference numbers come from",
        "",
        "* **57D torus bound:** we derive it from the published Schütt 56D input "
        "`5^28/(2^28*29^4)` using separated offsets and colored layers. "
        "The 57D number is our calculation, not a result attributed to Schütt. "
        "[Derivation](../code_routes/TORUS_LAYER_BENCHMARKS.md) and "
        "[exact radius arithmetic](propagated_benchmarks.py).",
        "* **Gamma72 sections:** we embed primitive laminated normals into the scaled "
        "Leech block of Nebe's parent and verify them in "
        "[verify_gamma72_leech_sections.py](verify_gamma72_leech_sections.py). "
        "For rank `k`, `delta^2=2^(72-k)/det(K)`. In 59D, `det(K)=2^23`, "
        "hence `delta=2^18=262144`. This specific 59D value was calculated here; "
        "it is not cited as a numerical result of Nebe's paper.",
        "* **Stronger 59D classical section:** Nebe's Section 4.2 supplies a primitive "
        "`sqrt(2)*Leech` section. The catalogue's classical `kappa13` embeds primitively "
        "in Leech with Gram determinant 972, as verified from both published coordinate matrices. "
        "The scaled normal and its 59D perpendicular have determinant `2^13*972`. "
        "The perpendicular has minimum exactly 8 because it contains Nebe's minimum-8 "
        "48D section. Hence `delta=2^22/(9*sqrt(3))=269064.727051…`, and our gain is "
        "**2.028457%**. This is a derivation from N12 and K13, not a density quoted "
        "from either source. [Proof and source survey](../../../constructions/dim59/README.md#literature-derived-comparison).",
        "* **65D depth-8 extension:** the existence proof using the Ne64 theta series "
        "was developed here: [NE64_THETA_COUNT.md](../code_routes/NE64_THETA_COUNT.md). "
        "Its center density is `3^(33/2)/2`. The explicit deeper coset gives "
        "`3^17/sqrt(10)`; the internal ratio is `sqrt(6/5)`.",
        "* **Additional propagated bounds:** products, torus colorings and averaged "
        "ball-window projections are retained in "
        "[propagated_benchmarks.json](propagated_benchmarks.json). Their numerical "
        "consequences are marked as our computations, even when their parent is published.",
        "",
        "## Source details",
        "",
    ]
    for sid, source in sources.items():
        lines += [
            f"**{sid}. {source['authors']}.** {link(source['title'], source['url'])}. "
            f"{source['publication']}. **Locator:** {source['locator']}.",
            "",
        ]
        if source.get("archived_source"):
            lines += [f"[Archived source]({source['archived_source']}).", ""]
        if source.get("limitation"):
            lines += [source["limitation"], ""]
    lines += [
        "## Replay and claim status",
        "",
        "```bash",
        "OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/propagated_benchmarks.py",
        "python constructions/dim59/verify.py",
        "python constructions/dim59/verify_independent.py",
        "python constructions/dim59/verify_external_comparison.py",
        "python constructions/dim59/verify_literature_section.py",
        "```",
        "",
        "The first command regenerates this page, the exact comparison arithmetic and "
        "source attribution from [benchmark_sources.json](benchmark_sources.json). "
        "It consumes previously verified construction data; it does not replace their "
        "minimum proofs. The next two commands replay the 59D construction and its "
        "internal reference independently. The fourth checks the historical reported-guarantee "
        "interval using exact rational powers; its source transcription is an input. "
        "The last verifies the classical kappa13 embedding and the stronger derived 59D comparison.",
        "",
        "**Status:** certified constructions with named published and derived comparisons. "
        "The dated literature review registers 59D as best-known. Other rows here establish "
        "the specified numerical comparisons. Sources and proved guarantees determine the "
        "comparison; no previous-record registry is required.",
        "",
    ]
    return "\n".join(lines)
