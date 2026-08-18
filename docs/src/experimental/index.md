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

There are two ways to set up a system: a text input file, or Julia code
calling the packing engine directly.

- [Input files](input_files.md): the global keywords and `structure ... end
  structure` block syntax.
- [Constraints](constraints.md): the geometric regions (`box`, `sphere`,
  `plane`, `cylinder`, `ellipsoid`, ...) that structures can be packed into
  — shared syntax between input files and the Julia API.
- [Defining systems in Julia](julia_api.md): building and packing a system
  with `structure_type`/`PackmolSystem`/`packmol`, without a text file.
- [Concentration Unit Conversion](@ref): part of `PackmolInputCreator`,
  which helps generate Packmol input files for common solution setups
  (solute + solvent, cosolvent, ions).

Everything else about the native packing engine is listed, undocumented
beyond docstrings, on the [Reference](reference.md) page. A proper usage
guide for the optimization loop itself will be added here as the engine
matures.
