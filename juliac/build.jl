#
# Builds the `packmol` executable with `juliac` (JuliaC.jl): a trimmed, relocatable
# bundle, with the executable in `<outdir>/bin` and the libraries it needs in
# `<outdir>/lib` (on Windows everything is in `<outdir>/bin`).
#
# Usage, from the root of the repository:
#
#     julia --project=juliac -e 'using Pkg; Pkg.instantiate()'
#     julia --project=juliac juliac/build.jl [outdir]
#
# `outdir` defaults to `juliac/build`. Requires Julia ≥ 1.13 (with 1.12 the bodies of
# spawned tasks are not compiled into the executable) and a C compiler (gcc or clang on
# Linux and macOS; on Windows JuliaC downloads MinGW). See the `juliac --trim`
# section of CLAUDE.md for what the code must satisfy to compile.
#
using JuliaC: ImageRecipe, LinkRecipe, BundleRecipe, compile_products, link_products, bundle_products

const repo = dirname(@__DIR__)
const outdir = abspath(isempty(ARGS) ? joinpath(@__DIR__, "build") : ARGS[1])

VERSION >= v"1.13" || error("building the executable requires Julia ≥ 1.13 (running $VERSION)")

image = ImageRecipe(
    output_type="--output-exe",
    file=joinpath(@__DIR__, "app.jl"),
    project=repo,
    trim_mode="safe",
    verbose=true,
)
link = LinkRecipe(
    image_recipe=image,
    outname=joinpath(outdir, "packmol"), # moved to bin/ (".exe" added on Windows)
    rpath="@bundle",
)
bundle = BundleRecipe(link_recipe=link, output_dir=outdir)

rm(outdir; force=true, recursive=true)
compile_products(image)
link_products(link)
bundle_products(bundle)
println("packmol executable built in $outdir")
