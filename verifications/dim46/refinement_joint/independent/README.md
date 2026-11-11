# Independent 46D border radial-refit audit

Read `AUDIT.md` for the complete proof review and scope. Run `python verify_all.py`
with Python 3.11+ and python-flint 0.9.0. Never use Python -O. No optimizer,
network, or files outside this package are needed.

All source packages were read only. Frozen inputs and all outputs are local to
this new independent package. The previous border_joint certificates remain
unchanged. `audit_result.json` records the aggregate independent outcome.
