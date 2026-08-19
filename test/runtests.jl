using TestItemRunner: @run_package_tests, @testitem

# Unitful defaults to printing exponents as Unicode superscripts (e.g. `mol⁻¹`)
# on macOS and as ASCII (`mol^-1`) everywhere else (see `Sys.isapple()` in
# Unitful/src/display.jl), unless this variable is set. Pin it so `show`-based
# tests (e.g. SolutionBoxUSC's) produce the same text on every CI platform.
ENV["UNITFUL_FANCY_EXPONENTS"] = "false"

@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(Packmol)
end

@testitem "Doctests" begin
    using Documenter: doctest
    doctest(Packmol)
end