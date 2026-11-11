# 21D Gaussian refinement release checks

The new registry value is0.004985977000. The full primary and independent proof chains are replayed from the staged package; the independent chain includes exact rational Sturm for both global rays, complete angular/cyclic and pair domains, source hashes, comparison arithmetic and five negative controls. See the complete replay log and audit report.

Historical proof artifacts are preserved byte for byte except origin_slack/verify_algebra.py, whose allowed registration paths now include the stronger Gaussian refinement. Its own sharp historical bound, exact normalization and Schur checks remain unchanged and were freshly replayed. Earlier historical audits are reused, not described as a new full historical aggregate. The earlier interrupted eight-dimension aggregate remains NOT PASSED as previously documented.

Changed integration tests cover registry/source hash rejection, exact proof subprocess commands, serial/parallel runner equivalence and failure propagation. Mocked subprocess tests are not mathematical replays. The comparison generator checks70/100-digit agreement and exact rational percentages/gap identities. All other registered rows are unchanged; the21D idea label and proof link now describe the Gaussian cyclic method.

Proof replay needs python-flint0.9.0, with assertions enabled. Discovery dependencies are separate. Published lower constructions and pinned upper references remain cited inputs, not independently reconstructed packings or exhaustive worldwide-priority claims.
