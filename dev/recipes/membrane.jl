#
# Images of the examples of the Membrane recipe (docs/src/recipes/membrane_system.md): the
# packed bilayers seen from the side, with the edges of the periodic cell, the lipids
# colored by type. Run from the root of the repository:
#
#     julia --project=docs -t 14 dev/recipes/membrane.jl
#
using Packmol
include(joinpath(@__DIR__, "..", "tutorials", "render.jl"))

const images = joinpath(@__DIR__, "..", "..", "docs", "src", "recipes", "images")
const work = mktempdir()

const view = (
    reps = [
        (selection = "water", style = "Licorice 0.2 12 12", color = "Name"),
        (selection = "not water", style = "Licorice 0.3 12 12", color = "ResName"),
    ],
    rotate = [(:x, -90), (:y, 20), (:x, 10)],
    cell = :origin,
)

for (name, system) in (
    "membrane_mixture" => Membrane(
        lipids = ["POPC", "POPE", "POPS", "POPG"], lipid_molar_ratio = [4, 3, 2, 1], total_lipids = 200,
    ),
    "membrane_dppc_dopc" => Membrane(lipids = ["DPPC", "DOPC"], total_lipids = 200),
)
    output = joinpath(work, "$name.pdb")
    t = @elapsed packed = packmol(system; output, seed = 1)
    println("$name: $(packed.status) in $(round(t; digits = 1)) s")
    render(output, joinpath(images, "$name.webp"); view...)
end
