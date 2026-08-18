import Pkg
Pkg.add("Documenter")
using Documenter
using Packmol
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [Packmol],
    sitename = "Packmol.jl",
    format = Documenter.HTML(top_menu = true),
    pages = [
        "Stable" => [
            "Home" => "index.md",
        ],
        "Experimental" => [
            "Home" => "experimental/index.md",
            "Input files" => "experimental/input_files.md",
            "Constraints" => "experimental/constraints.md",
            "Defining systems in Julia" => "experimental/julia_api.md",
            "Concentration units" => "concentration_units.md",
            "Reference" => "experimental/reference.md",
        ],
    ],
)
deploydocs(
    repo = "github.com/m3g/Packmol.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
