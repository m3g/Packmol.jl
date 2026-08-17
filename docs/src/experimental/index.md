# Packmol.jl (experimental)

!!! warning "Experimental"
    This section documents an in-development, pure-Julia reimplementation of
    Packmol that lives alongside the stable [`run_packmol`](@ref) wrapper
    described in the main documentation. The API and behavior described here
    can change at any time and are not yet covered by the same stability
    guarantees as the rest of the package.

## What this is

Packmol.jl is being rewritten natively in Julia, aiming for feature parity
with the classic Fortran [Packmol](http://github.com/m3g/packmol) package
while enabling new features, better maintainability, and closer integration
with the Julia ecosystem (for example, [CellListMap.jl](https://github.com/m3g/CellListMap.jl)
for efficient distance computations and [SPGBox.jl](https://github.com/m3g/SPGBox.jl)
as the optimizer).

This is not a wrapper around the Fortran binary: it reads the same kind of
input file, builds the packing problem, and solves it entirely in Julia.

## What's here today

There is no usage guide for the native packing engine yet; its public
functions are listed, undocumented beyond their docstrings, on the
[Reference](reference.md) page.

The one part of the new implementation with a narrative guide so far is
`PackmolInputCreator`, which helps generate Packmol input files for common
solution setups (solute + solvent, cosolvent, ions):

- [Concentration Unit Conversion](@ref)

More pages will be added here as the native packing engine matures.
