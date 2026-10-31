@testitem "compute_cm_bounds: extrapolated range of few molecules" setup=[TestFiles] begin
    # The range of the centers of n molecules underestimates the region they can
    # occupy: the bounds are widened by range/(n - 1) on each side, and clamped to
    # the cell under PBC.
    using Packmol
    using StaticArrays
    using LinearAlgebra

    water = TestFiles.test_file("structure_files/water.pdb")
    st = structure_type(water; number=6)
    sys = PackmolSystem([st]; output="")
    for (i, x) in enumerate(0.0:2.0:10.0) # range 10, n = 6: margin 2
        sys.molecule_positions[i] = Packmol.MoleculePosition(SVector(x, 0.0, -x), zero(SVector{3,Float64}))
    end
    lo, hi = Packmol.compute_cm_bounds(sys; extrapolate=false)
    @test lo[1] ≈ [0.0, 0.0, -10.0]
    @test hi[1] ≈ [10.0, 0.0, 0.0]
    lo, hi = Packmol.compute_cm_bounds(sys)
    @test lo[1] ≈ [-2.0, 0.0, -12.0]
    @test hi[1] ≈ [12.0, 0.0, 2.0]

    # Under PBC, the widened bounds do not go beyond the cell (here, from -11 to 11)
    sys = PackmolSystem([st]; output="", unitcell=Matrix{Float64}(Diagonal([22.0, 22.0, 22.0])))
    for (i, x) in enumerate(0.0:2.0:10.0)
        sys.molecule_positions[i] = Packmol.MoleculePosition(SVector(x, 0.0, -x), zero(SVector{3,Float64}))
    end
    lo, hi = Packmol.compute_cm_bounds(sys)
    @test lo[1] ≈ [-2.0, 0.0, -11.0]
    @test hi[1] ≈ [11.0, 0.0, 2.0]
end

@testitem "movebad!: one quota for the whole system" setup=[TestFiles] begin
    # Ten structure types of 6 molecules must not move more molecules than a
    # single type of 60 would (movefrac of the free molecules), even when every
    # molecule is bad.
    using Packmol
    using Random

    water = TestFiles.test_file("structure_files/water.pdb")
    sts = [structure_type(water; number=6, constraints=[InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])]) for _ in 1:10]
    sys = PackmolSystem(sts; output="", check=true, seed=1)
    Packmol.packmol(sys)
    fmol = ones(sys.nmols)
    free_mol_indices = collect(1:sys.nmols)
    mol_structure_type = Packmol._build_mol_structure_type(sys)
    rng = Random.Xoshiro(1)
    for _ in 1:20
        moved = Packmol.movebad!(sys, fmol, free_mol_indices, mol_structure_type, rng; movefrac=0.05)
        @test 1 <= length(moved) <= 3
        @test allunique(moved)
    end
end
