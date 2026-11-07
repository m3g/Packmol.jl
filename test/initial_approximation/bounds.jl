@testitem "initial approximation in the truncated octahedron" setup=[TestFiles] begin
    # Regression test: under PBC, the centers of mass of molecules that satisfy their
    # constraints can be stored at a periodic image, outside the bounds of their type.
    # The setup used to compute the bounds from those positions, and then moved the
    # molecules outside them to random points of the bounds, regardless of their
    # constraints, leaving molecules far from their regions in the initial approximation
    using Packmol
    using LinearAlgebra: norm
    water = TestFiles.test_file("structure_files/water.pdb")
    inner = structure_type(water; number = 100, constraints = [InsideSphere([0.0, 0.0, 0.0], 10.0)])
    outer = structure_type(water; number = 1200, constraints = [OutsideSphere([0.0, 0.0, 0.0], 13.0)])
    system = PackmolSystem([inner, outer];
        output = "", tolerance = 2.0, unitcell = octahedral_unitcell(Float64, 40.0), check_initial_placement = true, seed = 3,
    )
    redirect_stdout(devnull) do
        packmol(system)
    end
    atoms = get_atoms(system)
    r = [norm([a.x, a.y, a.z]) for a in atoms]
    @test maximum(r[1:300]) < 10.1
    @test minimum(r[301:end]) > 12.9
end

@testitem "_clamp_molecules_to_bounds! keeps molecules that satisfy the constraints" setup=[TestFiles] begin
    using Packmol
    using StaticArrays
    using Random: Xoshiro
    water = TestFiles.test_file("structure_files/water.pdb")
    st = structure_type(water; number = 2, constraints = [InsideSphere([0.0, 0.0, 0.0], 20.0)])
    system = PackmolSystem([st]; output = "")
    Packmol._center_reference_coordinates!(system)
    angles = zero(SVector{3,Float64})
    # Both outside the bounds: the first satisfies the constraint, the second does not
    system.molecule_positions[1] = Packmol.MoleculePosition(SVector(15.0, 0.0, 0.0), angles)
    system.molecule_positions[2] = Packmol.MoleculePosition(SVector(50.0, 0.0, 0.0), angles)
    lo, hi = [SVector(-5.0, -5.0, -5.0)], [SVector(5.0, 5.0, 5.0)]
    redirect_stdout(devnull) do
        Packmol._clamp_molecules_to_bounds!(system, [1, 2], lo, hi, Xoshiro(1))
    end
    @test system.molecule_positions[1].cm == SVector(15.0, 0.0, 0.0)
    @test all(lo[1] .<= system.molecule_positions[2].cm .<= hi[1])
end

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
