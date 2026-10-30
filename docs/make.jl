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
        "Tutorials" => [
            "Overview" => "tutorials/index.md",
            "Water and urea" => "tutorials/mixture.md",
            "Water–chloroform interface" => "tutorials/interface.md",
            "Lipid bilayer" => "tutorials/bilayer.md",
            "Spherical vesicle" => "tutorials/spherical.md",
            "Solvated protein" => "tutorials/solvprotein.md",
        ],
        "Experimental" => [
            "Home" => "experimental/index.md",
            "Input files" => "experimental/input_files.md",
            "Constraints" => "experimental/constraints.md",
            "Periodic boundary conditions" => "experimental/periodic_boundary_conditions.md",
            "Defining systems in Julia" => "experimental/julia_api.md",
            "Recipes" => "experimental/recipes.md",
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
    # Pull requests from branches of this repository are published to previews/PR<number>
    push_preview = true,
)
