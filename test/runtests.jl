using TestItemRunner: @run_package_tests, @testitem
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(Packmol)
end

@testitem "Doctests" begin
    using Documenter: doctest
    doctest(Packmol)
end