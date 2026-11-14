#!/usr/bin/env bash
set -euo pipefail
packing_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
packing_cxx="${CXX:-g++}"
mkdir -p "$packing_root/bin"
"$packing_cxx" -O3 -std=c++17 "$packing_root/contact_walk/shell.cpp" -o "$packing_root/bin/shell_discovery"
"$packing_cxx" -O3 -std=c++17 "$packing_root/mixed_codes/enumerate.cpp" -o "$packing_root/bin/mixed_enumerate"
"$packing_cxx" -O3 -std=c++17 "$packing_root/mixed_codes/search_parity.cpp" -o "$packing_root/bin/mixed_search"
"$packing_cxx" -O3 -std=c++17 "$packing_root/seed_deformations/periodic44_shell.cpp" -o "$packing_root/bin/periodic44_shell"
printf '%s\n' 'Built discovery and finite-code tools. The floating lattice tools do not certify empty searches.'
