"""Tests of the installed `packmol` command (run with pytest after installing the wheel)."""

import subprocess

WATER_XYZ = """\
3
water molecule
O   0.000   0.000   0.000
H   0.957   0.000   0.000
H  -0.240   0.927   0.000
"""

INPUT = """\
tolerance 2.0
filetype xyz
output output.xyz

structure water.xyz
  number 10
  inside box 0. 0. 0. 20. 20. 20.
end structure
"""


def _setup(tmp_path):
    (tmp_path / "water.xyz").write_text(WATER_XYZ)
    (tmp_path / "input.inp").write_text(INPUT)


def test_stdin(tmp_path):
    """packmol < input.inp, as with Fortran Packmol."""
    _setup(tmp_path)
    with open(tmp_path / "input.inp") as f:
        result = subprocess.run(["packmol"], stdin=f, capture_output=True, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    lines = (tmp_path / "output.xyz").read_text().splitlines()
    assert int(lines[0]) == 30


def test_input_option(tmp_path):
    """packmol -i input.inp, from another directory (paths relative to the input)."""
    _setup(tmp_path)
    result = subprocess.run(
        ["packmol", "-i", str(tmp_path / "input.inp")], capture_output=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "output.xyz").exists()


def test_version():
    import packmol

    result = subprocess.run(["packmol", "--version"], capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.strip() != ""
    assert packmol.__version__
