import Pkg
Pkg.add("Documenter")
using Documenter
using Packmol
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [Packmol],
    sitename = "Packmol.jl",
    pages = [
        "Home" => "index.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/Packmol.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
