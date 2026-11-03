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
