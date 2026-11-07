@testitem "_random_cm_in_cell: draws inside the periodic cell" begin
    # The draws in a box larger than the cell are those inside the cell: here the
    # bounding box of a truncated octahedron, of which 46% is outside it.
    using Packmol
    using StaticArrays
    using LinearAlgebra
    using Random

    uc = Packmol.octahedral_unitcell(Float64, 186.0)
    center = SVector(0.0, 0.0, 0.0)
    lo, hi = SVector(-93.0, -108.9, -113.2), SVector(93.0, 109.4, 113.5)
    rng = Xoshiro(1)
    cms = [Packmol._random_cm_in_cell(rng, lo, hi - lo, true, uc, center, true) for _ in 1:2000]
    @test all(cm -> Packmol._wrap_cm(cm, uc, center, true) == cm, cms)
    @test all(cm -> all(lo .<= cm .<= hi), cms)
    # Without PBC, the box
    cms = [Packmol._random_cm_in_cell(rng, lo, hi - lo, false, uc, center, true) for _ in 1:2000]
    @test all(cm -> all(lo .<= cm .<= hi), cms)
    @test count(cm -> Packmol._wrap_cm(cm, uc, center, true) != cm, cms) > 500
end

@testitem "_random_angles: uniform rotations" begin
    using Packmol
    using StaticArrays
    using Statistics: mean
    using Random: Xoshiro
    # The components of a rotated unit vector are uniform in [-1, 1] for uniform
    # rotations (mean absolute value 1/2). Uniform Euler angles, which this replaced,
    # gave 0.405 for the z component of the rotated z axis
    rng = Xoshiro(1)
    R = [Packmol.eulermat(Packmol._random_angles(rng, SVector{3,Float64})) for _ in 1:100_000]
    for axis in 1:3, component in 1:3
        @test abs(mean(abs(r[component, axis]) for r in R) - 0.5) < 0.01
    end
    # In two dimensions, the angles are uniform
    a = [Packmol._random_angles(rng, SVector{2,Float64}) for _ in 1:1000]
    @test all(all(0 .<= x .< 2π) for x in a)
end
