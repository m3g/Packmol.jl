# Polymer in DMF

A chain of polyacrylamide solvated by N,N-dimethylformamide (DMF), in a periodic box: a
system of the [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl)
[examples](https://github.com/m3g/ComplexMixturesExamples), where the solvation of the
polymer by DMF is studied.

Files: [`polyacrylamide.pdb`](files/polymer_dmf/polyacrylamide.pdb),
[`dmf.pdb`](files/polymer_dmf/dmf.pdb), and the input file
[`polymer_dmf.inp`](files/polymer_dmf/polymer_dmf.inp).

The chain has five monomers (55 atoms). In Packmol, molecules are rigid: the chain keeps
the conformation of its file, and only its position and orientation change. For long
polymers, the conformation must come from elsewhere (an equilibrated simulation, or a
tool that builds chains), or the chain can be packed in parts, as different molecules
held together by constraints.

## The number of molecules

The box is a cube with sides of 40 Å. DMF (73.09 g/mol) has a density of 0.944 g/mL, that
is, 0.944 × 0.6022/73.09 = 0.00778 molecules per Å³. The chain occupies about 400 Å³, so
the box holds (40³ − 400) × 0.00778 ≈ 495 molecules of DMF.

## Command line

The input file:

```
#
# A polyacrylamide chain in N,N-dimethylformamide (DMF, 0.944 g/mL), in a
# periodic cubic box with sides of 40 Å.
#
tolerance 2.0
output polymer_dmf.pdb
pbc -20. -20. -20. 20. 20. 20.

structure polyacrylamide.pdb
  number 1
end structure

structure dmf.pdb
  number 495
end structure
```

The periodic box confines all the molecules: no constraint is needed. Run it with:

```
packmol polymer_dmf.inp
```

The packing takes a fraction of a second.

## Julia

```@setup polymer
for file in ("polyacrylamide.pdb", "dmf.pdb")
    cp(joinpath("files", "polymer_dmf", file), file; force = true)
end
```

```@example polymer
using Packmol, PDBTools

polymer = structure_type("polyacrylamide.pdb"; number = 1)
dmf = structure_type("dmf.pdb"; number = 495)
system = PackmolSystem([polymer, dmf];
    output = "polymer_dmf.pdb",
    tolerance = 2.0,
    unitcell = [40. 0. 0.; 0. 40. 0.; 0. 0. 40.], # centered at the origin
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

The packed system, with the `visualize` function of
[PDBTools.jl](https://github.com/m3g/PDBTools.jl): the polymer, and DMF as lines, with the
edges of the periodic box (drag to rotate, scroll to zoom):

```@example polymer
visualize(get_atoms(system),
    "not resname DMF" => (style = :ballandstick,),
    "resname DMF" => (style = :lines,);
    unitcell = UnitCell(system.unitcell),
    unitcell_origin = -sum(eachcol(system.unitcell)) / 2,
)
```
