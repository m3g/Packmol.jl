```@meta
CollapsedDocStrings = true
```

# Reference

!!! warning "Experimental"
    These are the public functions of the native Julia packing engine.
    Signatures and behavior may still change. Recipes' functions are
    documented on their own [Recipes](recipes.md) page.

## Packing engine

```@docs
packmol
structure_type
PackmolSystem
get_atoms
write_output
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
