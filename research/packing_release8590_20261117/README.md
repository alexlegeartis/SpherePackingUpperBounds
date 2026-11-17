# Six audited packing improvements, 85–90 dimensions

Method: **Correlated Mordell-Weil dual-normal sections with optimized simplex translate clusters**.

These are periodic packing lower bounds, not lattice-only assertions. The exact source-to-period construction and a separate Smith/coarea audit are included. See [proof and attribution](PROOF.md) and [machine-readable results](RESULTS.json).

Install `requirements.txt`, then run:

```
python run_all.py
```

Do not use `python -O`. The replay performs all six primary checks and the three independent audit batches. It needs no search solver or network. The inherited parent-minimum theorem and the shared source-reconstruction routine are stated explicitly in the proof.

The display floors and comparison gains are taken from the audited results. Baseline comparisons are pinned; they do not establish exhaustive worldwide priority.
