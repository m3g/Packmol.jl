```@meta
CollapsedDocStrings = true
```

# Reference

These are the public functions of Packmol.jl. The functions of the recipes are
documented on the [Recipes](../recipes/index.md) page, and the volume estimates,
`guess_volume` and `free_volume`, on the [Volume estimates](@ref) page.

## Packing engine

```@docs
packmol
structure_type
PackmolSystem
get_atoms
write_output
```

## Placement on a grid

A direct port to Julia of the placement algorithm of
[Bentopy](https://github.com/marrink-lab/bentopy): see [Placement on a grid (Bentopy)](@ref).

```@docs
voxel_place
```

## Tools for building systems

Public, but not exported: call them as `Packmol.cross_section`, `Packmol.builtin_lipid`.

```@docs
Packmol.cross_section
Packmol.builtin_lipid
```

## Shape of molecules

```@docs
shape_descriptors
select_conformations
```

## Box constraints

```@docs
InsideBox
OutsideBox
InsideCube
OutsideCube
```

## Dodecahedral periodic box

```@docs
dodecahedral_unitcell
triclinic_to_dodecahedral
dodecahedral_to_triclinic
```

## Octahedral periodic box

```@docs
octahedral_unitcell
triclinic_to_octahedral
octahedral_to_triclinic
```

## Legacy (Fortran) Packmol and the file browser

```@docs
run_packmol_legacy
Packmol.pick_file
```
