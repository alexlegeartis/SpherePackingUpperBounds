# Independent 21D Gaussian cyclic audit

PASS: unrestricted sphere-packing density Delta_21 < 0.004985977000.

Use Python 3.12 and python-flint==0.9.0, with assertions enabled:

    python -m pip install -r requirements.txt
    python verify_all.py

The self-contained replay needs no optimizer, network, source producer directory, or previous radial witness. It checks frozen-input hashes, dimension-specific angular and cyclic domains, exact rational Fourier/Sturm signs, PSD origin normalization, interval pair coupling, rational outward rounding, comparison arithmetic, and negative controls. Exact Sturm calculations can take several minutes.

AUDIT.md explains the analytic packing theorem and every limitation. full_audit.log is the final full-chain evidence; audit_exit_code.txt is 0. Both angular_cyclic_audit.json and audit.json are required. The inherited_dependency text in the radial-only result is discharged by the completed angular stage. All source binding evidence and comparison provenance are included.

SHA-256 checksums of this delivered audit are in AUDIT_SHA256SUMS. Initial preparatory logs are retained with initial_ prefixes and do not supersede the clean final run.
