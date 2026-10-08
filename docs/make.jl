using Documenter
using Packmol
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [Packmol],
    sitename = "Packmol.jl",
    format = Documenter.HTML(
        top_menu = true,
        # The interactive 3D views of the tutorials are embedded in the page
        size_threshold_ignore = ["experimental/tutorials.md"],
    ),
    pages = [
        "Stable" => [
            "Home" => "index.md",
        ],
        "Experimental" => [
            "Home" => "experimental/index.md",
            "Tutorials" => "experimental/tutorials.md",
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
# Remove the structures packed by the tutorials, which are not part of the docs
for file in readdir(joinpath(@__DIR__, "build", "experimental"); join = true)
    endswith(file, ".pdb") && rm(file)
end
deploydocs(
    repo = "github.com/m3g/Packmol.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
