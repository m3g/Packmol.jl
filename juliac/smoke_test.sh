#!/bin/bash
#
# Packs 10 water molecules with the given executable (default juliac/build/bin/packmol):
# a check that the bundle runs on a system, e.g. in a container of an old distribution:
#
#     docker run --rm -v "$PWD":/work -w /work rockylinux:8 bash juliac/smoke_test.sh
#
set -euo pipefail
exe=$(cd "$(dirname "${1:-juliac/build/bin/packmol}")" && pwd)/$(basename "${1:-juliac/build/bin/packmol}")
dir=$(mktemp -d)
cat > "$dir/water.xyz" <<'XYZ'
3
water
O   0.000   0.000   0.000
H   0.957   0.000   0.000
H  -0.240   0.927   0.000
XYZ
cat > "$dir/input.inp" <<'INP'
tolerance 2.0
output output.xyz
structure water.xyz
  number 10
  inside box 0. 0. 0. 20. 20. 20.
end structure
INP
(cd "$dir" && "$exe" < input.inp > packmol.log) || { cat "$dir/packmol.log"; exit 1; }
natoms=$(head -1 "$dir/output.xyz" | tr -d ' ')
[ "$natoms" = 30 ] || { echo "unexpected output: $natoms atoms"; exit 1; }
echo "smoke test passed ($(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME" || uname -s))"
rm -rf "$dir"
