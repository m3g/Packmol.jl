"""Hatchling build hook: puts the prebuilt `packmol` bundle into the wheel.

The bundle (the directory built by juliac/build.jl, with `bin/`, `lib/` and `share/`) is
taken from the environment variable PACKMOL_BUNDLE and copied to `packmol/binaries/`.

Wheels cannot hold symbolic links (they would be stored as copies, and the libraries of
the bundle have up to three names each), so the links are left out and listed in
`packmol/binaries/symlinks.txt`; `packmol.cli` creates them when first run.

The platform tag is PACKMOL_WHEEL_PLATFORM if set, otherwise:
- Linux: `manylinux_2_28_<arch>` (the bundle is built in a manylinux_2_28 container, and
  the workflow checks that it needs no newer glibc);
- macOS: `macosx_<major>_<minor>_<arch>`, from the minimum macOS version of its files;
- Windows: `win_amd64`.
"""

import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


def _macos_min_version(bundle: Path) -> tuple[int, int]:
    """Largest minimum macOS version (LC_BUILD_VERSION minos) of the Mach-O files."""
    vmax = (11, 0)
    for f in bundle.rglob("*"):
        if f.is_symlink() or not f.is_file():
            continue
        out = subprocess.run(["otool", "-l", str(f)], capture_output=True, text=True)
        for m in re.finditer(r"minos (\d+)\.(\d+)", out.stdout):
            v = (int(m.group(1)), int(m.group(2)))
            if v > vmax:
                print(f"minimum macOS {v[0]}.{v[1]}: {f.relative_to(bundle)}")
                vmax = v
    # From macOS 11 on, the minor version of the tag is 0
    return (vmax[0], 0) if vmax[0] >= 11 else vmax


def _platform_tag(bundle: Path) -> str:
    tag = os.environ.get("PACKMOL_WHEEL_PLATFORM")
    if tag:
        return tag
    machine = platform.machine().lower()
    if sys.platform.startswith("linux"):
        arch = {"amd64": "x86_64", "arm64": "aarch64"}.get(machine, machine)
        return f"manylinux_2_28_{arch}"
    if sys.platform == "darwin":
        major, minor = _macos_min_version(bundle)
        arch = "arm64" if machine in ("arm64", "aarch64") else "x86_64"
        return f"macosx_{major}_{minor}_{arch}"
    if sys.platform == "win32":
        return "win_amd64"
    raise RuntimeError(f"unsupported platform: {sys.platform} {machine}")


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        if self.target_name != "wheel":
            raise RuntimeError("only wheels can be built (there is no source distribution)")
        bundle_env = os.environ.get("PACKMOL_BUNDLE")
        if not bundle_env:
            raise RuntimeError("set PACKMOL_BUNDLE to the bundle built by juliac/build.jl")
        bundle = Path(bundle_env)
        if not bundle.is_absolute():
            # Relative to the root of the repository
            bundle = Path(self.root).parent / bundle
        exe = "packmol.exe" if sys.platform == "win32" else "packmol"
        if not (bundle / "bin" / exe).is_file():
            raise RuntimeError(f"executable not found: {bundle / 'bin' / exe}")

        target = Path(self.root) / "packmol" / "binaries"
        shutil.rmtree(target, ignore_errors=True)
        target.mkdir(parents=True)
        links = []
        for f in sorted(bundle.rglob("*")):
            rel = f.relative_to(bundle)
            if f.is_symlink():
                link = os.readlink(f)
                if os.path.isabs(link) or ".." in Path(link).parts:
                    raise RuntimeError(f"link out of its directory in the bundle: {rel} -> {link}")
                links.append(f"{rel.as_posix()} {link}")
            elif f.is_dir():
                (target / rel).mkdir(exist_ok=True)
            else:
                shutil.copy2(f, target / rel)
        (target / "symlinks.txt").write_text("".join(line + "\n" for line in links))

        build_data["pure_python"] = False
        build_data["tag"] = f"py3-none-{_platform_tag(bundle)}"
        # The bundle is not a Python package: include it as artifacts (it is gitignored)
        build_data["artifacts"].append("packmol/binaries/")

    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        shutil.rmtree(Path(self.root) / "packmol" / "binaries", ignore_errors=True)
