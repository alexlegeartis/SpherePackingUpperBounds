# Frozen midrange release evidence

This release registers dimensions25,28,32,34,36,38,39,40. Each has a full primary proof, a separately implemented completed audit and portable discovery/export code. The certificate inventory gives exact outward bounds and citations.

The attempted fresh aggregate was interrupted and **did not pass**. Previously completed full proofs are reused only through the exact file bindings in AUDIT_EQUIVALENCE.json; old proof subtrees and the completed52fa554 regression are identified in HISTORICAL_BINDING.json. All changed integration, registration, export and runner checks passed. RELEASE_CHECKS.json separates these scopes explicitly. Requirement files and the documented python-flint0.9.0 version match; original runtime binary hashes were not recorded and are not asserted.

Run `python comparison.py` at repository root for table checks. Run `python research/gaussian_midrange_release/verify_bindings.py` for the packaged evidence and integration checks. This is not a fresh mathematical replay. To replay every proof from scratch, run `python run_all.py` or `python run_all.py --jobs 4`. Different dimensions can run concurrently, while each dimension remains sequential and shared final checks wait for all chains. Numerical search outputs are never used as proof.
