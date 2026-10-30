#
# Images of the water drop tutorial (docs/src/tutorials/water_drop.md): packs the inputs
# of docs/src/tutorials/files/water_drop with the executable (juliac/build.jl), with and
# without `check`, and renders them.
#
#     julia dev/tutorials/water_drop.jl
#
include("render.jl")
const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
cp(joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "water_drop"), work; force=true)

function pack(name; check=false)
    input = read(joinpath(work, "$name.inp"), String)
    out = check ? "$(name)_check" : name
    check && (input = replace(input, "output $name.pdb" => "output $out.pdb\ncheck"))
    write(joinpath(work, "$out.inp"), input)
    run(pipeline(`$exe -i $(joinpath(work, "$out.inp"))`; stdout=devnull))
    return joinpath(work, "$out.pdb")
end

reps = [
    (selection = "resname GRA", style = "VDW 0.7 12", color = "ColorID 2"),
    (selection = "water", style = "Licorice 0.25 12 12"),
]
view = (rotate = [(:x, -75), (:y, 20)],)
render(pack("drop"; check=true), joinpath(images, "drop_check.webp"); reps, view...)
render(pack("drop"), joinpath(images, "drop.webp"); reps, view...)
render(pack("gaussian_drop"), joinpath(images, "gaussian_drop.webp"); reps, view...)
