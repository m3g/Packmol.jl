# Legacy (Fortran) Packmol

The legacy Packmol, the Fortran program, is still available from
[m3g/packmol](http://github.com/m3g/packmol) (see
[Migrating from legacy Packmol](@ref) for the differences). It can also be run from Julia, with
[`run_packmol_legacy`](@ref), which works as `packmol`. It requires the `Packmol_jll`
package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
```
