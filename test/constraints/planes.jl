@testitem "Plane constructors" begin
    @test AbovePlane([0, 0, 1], 5.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test AbovePlane(normal = [0, 0, 1], d = 5.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test AbovePlane(normal = [0, 0, 1], d = 5.0, weight = 2.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 2.0)
    @test BelowPlane([0, 0, 1], 5.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test BelowPlane(normal = [0, 0, 1], d = 5.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test BelowPlane(normal = [0, 0, 1], d = 5.0, weight = 2.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 2.0)
end

@testitem "Plane gradients" begin
    using ForwardDiff
    using StaticArrays
    # Test above plane: n = (0,0,1), d = 5 -> atom must have z >= 5
    # Point below the plane (should have penalty)
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    c = AbovePlane([0, 0, 1], 5.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Point above the plane (no penalty)
    x = SVector{3,Float64}(1.0, 2.0, 7.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Test below plane: n = (0,0,1), d = 5 -> atom must have z <= 5
    # Point above the plane (should have penalty)
    x = SVector{3,Float64}(1.0, 2.0, 7.0)
    c = BelowPlane([0, 0, 1], 5.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Point below the plane (no penalty)
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Test with tilted plane normal
    x = SVector{3,Float64}(1.5, 1.0, 0.5)
    c = AbovePlane([1, 1, 1] ./ sqrt(3.0), 2.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    c = BelowPlane([1, 1, 1] ./ sqrt(3.0), 0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
end
