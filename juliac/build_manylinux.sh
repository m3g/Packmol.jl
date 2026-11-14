#!/bin/bash
#
# Builds the bundle (juliac/build.jl) in a manylinux_2_28 container (AlmaLinux 8, glibc
# 2.28), so that the Linux executable runs on RHEL/Rocky/Alma 8 and newer, Ubuntu ≥ 20.04
# and Debian ≥ 11. Built on a newer system, it would need the glibc of that system.
#
# Run from the root of the repository, in the container (the CI workflow does this):
#
#     docker run --rm -v "$PWD":/work -w /work quay.io/pypa/manylinux_2_28_x86_64 \
#         bash juliac/build_manylinux.sh
#
# (quay.io/pypa/manylinux_2_28_aarch64 on arm64). JULIA_VERSION selects the Julia version.
#
set -euo pipefail
JULIA_VERSION=${JULIA_VERSION:-1.13.2}
arch=$(uname -m)
case $arch in
    x86_64) dir=x64 ;;
    aarch64) dir=aarch64 ;;
    *) echo "unsupported architecture: $arch" >&2; exit 1 ;;
esac
julia_dir=/tmp/julia-$JULIA_VERSION
if [ ! -x "$julia_dir/bin/julia" ]; then
    curl -fsSL "https://julialang-s3.julialang.org/bin/linux/$dir/${JULIA_VERSION%.*}/julia-$JULIA_VERSION-linux-$arch.tar.gz" |
        tar xz -C /tmp
fi
export PATH=$julia_dir/bin:$PATH
gcc --version | head -1
julia --project=juliac -e 'using Pkg; Pkg.instantiate()'
julia --project=juliac juliac/build.jl
