# Entry point of the `packmol` executable, compiled with `juliac` (see build.jl)
using Packmol
(@main)(ARGS) = Packmol.CLI.main(ARGS)
