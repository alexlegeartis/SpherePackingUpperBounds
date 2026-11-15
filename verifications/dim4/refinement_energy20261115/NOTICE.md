# Attribution and license

This extension builds on **Felpix-Studios/sphere-packing**,
<https://github.com/Felpix-Studios/sphere-packing>, commit
`bb97550ba47ad92bacfee567843888be16895a23`, particularly its preprint
*An improved upper bound for sphere packing density in four dimensions*.

The 277 files under `vendor/D4/` are unchanged upstream copies. Their
per-file provenance and SHA-256 hashes are in `vendor/UPSTREAM.json`.
The upstream Apache License 2.0 is retained in
[`vendor/LICENSE.Apache-2.0`](vendor/LICENSE.Apache-2.0). All upstream
copyright, attribution and license notices remain intact. Upstream did
not ship a separate NOTICE file at the pinned commit.

The upstream authors supplied the underlying scalar/kernel certificate,
count and graph refinements, original far-pair argument, interval engine
and complete base replay. The new contribution is the degree-16 angular
energy polynomial, extraction of unused angular Gram credits, its exact
and interval certificates, and the resulting stricter bound.

`boxes.py` adapts the upstream subdivision format and traversal. Its file
header identifies the sources and modifications. It is distributed under
Apache-2.0 as well. The other new files in this package are also made
available under Apache-2.0. This license statement does not relicense
unrelated files elsewhere in this repository.
