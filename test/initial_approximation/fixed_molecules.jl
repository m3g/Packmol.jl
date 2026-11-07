@testitem "overlaps_fixed: periodic images and avoid_overlap" setup=[TestFiles] begin
    # A trial molecule near a face of the periodic cell overlaps a fixed molecule near
    # the opposite face, through its periodic image. And the relocations of the packing
    # loop avoid the fixed molecules even with `avoid_overlap no`.
    using Packmol
    using StaticArrays
    using LinearAlgebra

    water = TestFiles.test_file("structure_files/water.pdb")
    fixed = structure_type(water; number=1, fixed=([10.5, 0.0, 0.0], [0.0, 0.0, 0.0]), center=:geometric)
    free = structure_type(water; number=1)
    sys = PackmolSystem([fixed, free]; output="", unitcell=Matrix{Float64}(Diagonal([22.0, 22.0, 22.0])))
    tol = sys.tolerance
    # (the reference coordinates are centered only when the packing starts)
    ref = sys.structure_types[2].reference_coordinates
    ref = ref .- Ref(sum(ref) / length(ref))
    buffer = similar(ref)
    angles = zero(SVector{3,Float64})
    fixed_sys, lo, hi = Packmol._build_overlap_check_system(sys, tol)
    @test !isnothing(fixed_sys)
    # 1 Å from the fixed molecule through the face at x = 11 (the cell is centered at
    # the origin), 21 Å from it in the cell
    @test Packmol.overlaps_fixed(Packmol.MoleculePosition(SVector(-10.5, 0.0, 0.0), angles), ref, fixed_sys, buffer, lo, hi, tol)
    @test Packmol.overlaps_fixed(Packmol.MoleculePosition(SVector(9.5, 0.0, 0.0), angles), ref, fixed_sys, buffer, lo, hi, tol)
    @test !Packmol.overlaps_fixed(Packmol.MoleculePosition(SVector(0.0, 0.0, 0.0), angles), ref, fixed_sys, buffer, lo, hi, tol)

    # Without PBC, the far molecule does not overlap
    sys = PackmolSystem([fixed, free]; output="")
    fixed_sys, lo, hi = Packmol._build_overlap_check_system(sys, tol)
    @test !Packmol.overlaps_fixed(Packmol.MoleculePosition(SVector(-10.5, 0.0, 0.0), angles), ref, fixed_sys, buffer, lo, hi, tol)
    @test Packmol.overlaps_fixed(Packmol.MoleculePosition(SVector(9.5, 0.0, 0.0), angles), ref, fixed_sys, buffer, lo, hi, tol)

    # avoid_overlap = false disables the check of the initial approximation only
    sys = PackmolSystem([fixed, free]; output="", avoid_overlap=false)
    @test isnothing(first(Packmol._build_overlap_check_system(sys, tol)))
    @test !isnothing(first(Packmol._build_overlap_check_system(sys, tol; avoid_overlap=true)))
end
