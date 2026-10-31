```@meta
CollapsedDocStrings = true
```

# Reference

These are the public functions of Packmol.jl. The functions of the recipes are
documented on the [Recipes](../recipes/index.md) page.

## Packing engine

```@docs
packmol
structure_type
PackmolSystem
get_atoms
write_output
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
