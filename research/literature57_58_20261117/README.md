# 57D/58D numerical literature comparator extraction

Checked 2026-11-17. This is a source-data extraction, **not a rigorous certification of the old LP benchmark and not a complete global-priority search**.

Source: Afkhami-Jeddi, Cohn, Hartman, de Laat, Tajdini, *High-dimensional sphere packing and the modular bootstrap* (2020), https://arxiv.org/abs/2006.02560 . Official dataset: https://hdl.handle.net/1721.1/125646 .

Downloaded `dataminusone.txt` from https://dspace.mit.edu/server/api/core/bitstreams/8651d342-73d4-4b91-bba8-7074de435e4a/content . Size 283480522 bytes; MD5 `2c062199a2305ae17f512dc00f68d735` matches the repository metadata. Dataset README identifies this as numerically optimized scaling dimensions and degeneracies with total central charge equal to dimension. Full source retained outside this small extraction package.

Conversion: the first scaling dimension is Delta1=r²/2 (paper Section 2.2); Eq. (2.18) gives density pi^(d/2)*(r/2)^d/Gamma(d/2+1). Therefore numerical density=(pi*Delta1/2)^(d/2)/Gamma(d/2+1). The script uses 80-decimal arithmetic, which does not turn numerical source values into rigorously bounded intervals.

Highest available truncations for these dimensions:

|d|N|Numerical LP comparator|Published repository upper|Percent below numerical comparator|
|--:|--:|--:|--:|--:|
|57|493|1.7231960044411017e-8|1.5718326e-8|8.783876%|
|58|496|1.1947773131075187e-8|1.073633e-8|10.139489%|

Paper Section 2.4 describes numerical solutions and the additional positivity check needed for rigor. Section 2.3 discusses conjectural lifting of the optimal eigenfunction into the Cohn–Elkies formulation. Do not label extracted source values as independently certified prior upper bounds. Cohn's public density table omits dimensions 57 and 58, so these are sourced from the full modular-bootstrap dataset, not that table.

Replay: `python extract.py /path/to/dataminusone.txt` (requires mpmath). `results.json` contains exact source decimals and computed comparison values.

Independent second extraction verified the file checksum, selected highest-N records, normalization and percentages. Source gap lines are 1595340 (57D) and 1597428 (58D). This audit validates the comparison calculation, not the underlying external optimizer.
