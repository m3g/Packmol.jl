using TestItemRunner: @run_package_tests, @testitem

# Unitful defaults to printing exponents as Unicode superscripts (e.g. `mol⁻¹`)
# on macOS and as ASCII (`mol^-1`) everywhere else (see `Sys.isapple()` in
# Unitful/src/display.jl), unless this variable is set. Pin it so `show`-based
# tests (e.g. SolutionBoxUSC's) produce the same text on every CI platform.
ENV["UNITFUL_FANCY_EXPONENTS"] = "false"

@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    # The Downgrade workflow's julia-downgrade-compat action copies the test-only
    # dependencies into [deps] (to keep their lower bounds in the Pkg.test sandbox),
    # which Aqua would report as stale.
    Aqua.test_all(Packmol; stale_deps = !haskey(ENV, "PACKMOL_DOWNGRADE_TEST"))
end
