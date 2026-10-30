# Restart files

`restart_from`/`restart_to`, on both `PackmolSystem` (whole system) and
`structure_type` (one structure type only), work as in the input file's
[Restarting a run](../input_files/restarting_a_run.md) section — a coordinate
file (PDB, XYZ or mmCIF, by extension) skips or shortcuts the initial placement
step:

```julia
sys = PackmolSystem(structure_types; output="box.pdb", tolerance=2.0, restart_to="box_restart.pdb")
packmol(sys)  # ... later, resume packing (or extend a run that stopped early):
sys2 = PackmolSystem(structure_types; output="box2.pdb", tolerance=2.0, restart_from="box_restart.pdb")
packmol(sys2)
```

The Julia API additionally accepts `restart_from` as an already-loaded
`Vector{<:Atom}` (from `PDBTools.read_pdb` (PDBTools ≥ 4), or built up any other way)
instead of a file path, skipping the file entirely:

```julia
using PDBTools: read_pdb
atoms = read_pdb("some_configuration.pdb")
sys3 = PackmolSystem(structure_types; output="box3.pdb", tolerance=2.0, restart_from=atoms)
```
