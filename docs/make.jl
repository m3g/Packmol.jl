using Documenter, DocumenterCodeBlocks
using Packmol
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [Packmol],
    sitename = "Packmol.jl",
    format = Documenter.HTML(top_menu = true),
    plugins = [CodeBlocks()],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Overview" => "tutorials/index.md",
            "Water and urea" => "tutorials/mixture.md",
            "Water–chloroform interface" => "tutorials/interface.md",
            "Lipid bilayer" => "tutorials/bilayer.md",
            "Spherical vesicle" => "tutorials/spherical.md",
            "Solvated protein" => "tutorials/solvprotein.md",
        ],
        "Input files" => "input_files.md",
        "Constraints" => "constraints.md",
        "Periodic boundary conditions" => "periodic_boundary_conditions.md",
        "Julia API" => "julia_api.md",
        "Recipes" => "recipes.md",
        "Concentration units" => "concentration_units.md",
        "Reference" => "reference.md",
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
