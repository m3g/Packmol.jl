"""Command-line entry point: runs the bundled `packmol` executable.

The wheel holds the bundle built by juliac (the executable and the Julia libraries) without
its symbolic links, which wheels cannot store; they are listed in `binaries/symlinks.txt`
and created here when first run. If the installation directory is not writable, the bundle
is copied to a cache directory of the user, and the links are created there.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_MARKER = ".symlinks-done"


def _make_links(bundle: Path) -> None:
    listing = bundle / "symlinks.txt"
    lines = listing.read_text().splitlines() if listing.exists() else []
    for line in lines:
        name, target = line.split(" ", 1)
        link = bundle / name
        if link.exists() or link.is_symlink():
            continue
        try:
            os.symlink(target, link)
        except FileExistsError:
            pass  # created by another process
        except OSError:
            if os.name != "nt":
                raise
            # Windows without the privilege of creating links: a copy
            shutil.copy2(link.parent / target, link)
    (bundle / _MARKER).touch()


def _cache_dir() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    try:
        from importlib.metadata import version

        v = version("packmol")
    except Exception:
        v = "unknown"
    return base / "packmol" / v


def get_bundle() -> Path:
    """The directory of the bundle, with its links in place."""
    bundle = Path(__file__).parent / "binaries"
    if not bundle.is_dir():
        raise FileNotFoundError(f"packmol bundle not found in {bundle}: broken installation")
    if (bundle / _MARKER).exists():
        return bundle
    try:
        _make_links(bundle)
        return bundle
    except OSError:
        pass
    # The installation is not writable: a copy in the cache of the user
    cached = _cache_dir()
    if (cached / _MARKER).exists():
        return cached
    cached.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(dir=cached.parent))
    try:
        shutil.copytree(bundle, tmp / "bundle", symlinks=True)
        # The copy keeps the permissions of the installation: make its directories writable
        for d, _, _ in os.walk(tmp / "bundle"):
            os.chmod(d, os.stat(d).st_mode | 0o700)
        _make_links(tmp / "bundle")
        try:
            os.replace(tmp / "bundle", cached)
        except OSError:
            if not (cached / _MARKER).exists():  # not made by another process meanwhile
                raise
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return cached


def get_binary_path() -> Path:
    """Path of the packmol executable."""
    exe = "packmol.exe" if os.name == "nt" else "packmol"
    path = get_bundle() / "bin" / exe
    if not path.is_file():
        raise FileNotFoundError(f"packmol executable not found: {path}")
    return path


def main() -> int:
    try:
        binary = get_binary_path()
    except (OSError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    args = [str(binary)] + sys.argv[1:]
    if os.name == "nt":
        return subprocess.run(args, check=False).returncode
    os.execv(args[0], args)
    return 1  # not reached


if __name__ == "__main__":
    sys.exit(main())
