# Legacy (Fortran) Packmol

The classical Fortran program can still be run from Julia, with
[`run_packmol_legacy`](@ref), which works as `packmol`. It requires the `Packmol_jll`
package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
```
