@testitem "octahedral_unitcell" begin
    uc = Packmol.octahedral_unitcell(Float64, 100.0)
    a, b, c, α, β, γ = Packmol._unitcell_abc_angles(uc)
    @test a ≈ 100.0
    @test b ≈ 100.0
    @test c ≈ 100.0
    @test α ≈ acosd(1 / 3)
    @test β ≈ acosd(-1 / 3)
    @test γ ≈ acosd(1 / 3)
    # matches the well-known GROMACS truncated-octahedron volume ratio:
    # ~76.98% of a cube's volume for the same minimum-image distance `d`
    # (vs. ~70.71% for the rhombic dodecahedron).
    @test abs(Packmol.det(uc)) / 100.0^3 ≈ 0.7698 atol = 1e-4
end

@testitem "parse_pbc_octahedral" begin
    using StaticArrays: SVector
    uc, center = Packmol.parse_pbc_octahedral(Float64, ["1.0", "2.0", "3.0", "100.0"], 3)
    @test uc ≈ Packmol.octahedral_unitcell(Float64, 100.0)
    @test center ≈ SVector(1.0, 2.0, 3.0)
    @test_throws ArgumentError Packmol.parse_pbc_octahedral(Float64, ["1.0", "2.0", "3.0", "100.0"], 2)
    @test_throws ArgumentError Packmol.parse_pbc_octahedral(Float64, ["1.0", "2.0", "3.0"], 3)
end

@testitem "pbc octahedral input keyword" begin
    using Packmol: read_packmol_input
    dir = dirname(Packmol.src_dir * "/../test/input_files/water_box.inp")
    original = read(Packmol.src_dir * "/../test/input_files/water_box.inp", String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))
        file = joinpath(tmp, "octa.inp")
        write(file, original * "\npbc octahedral 0. 0. 0. 100.\n")
        sys = read_packmol_input(file)
        @test !isnothing(sys.unitcell)
        a, b, c, α, β, γ = Packmol._unitcell_abc_angles(sys.unitcell)
        @test a ≈ b ≈ c ≈ 100.0
        @test α ≈ γ ≈ acosd(1 / 3)
        @test β ≈ acosd(-1 / 3)
        @test sys.unitcell_center == zeros(3)
    end
end

@testitem "triclinic_to_octahedral / octahedral_to_triclinic round trip" begin
    using StaticArrays: SVector
    using LinearAlgebra: norm
    using Random

    d = 100.0
    uc = Packmol.octahedral_unitcell(Float64, d)
    center = SVector(0.0, 0.0, 0.0)

    Random.seed!(1)
    for _ in 1:2000
        x = SVector(d * (rand() - 0.5), d * (rand() - 0.5), d * (rand() - 0.5))
        xt = Packmol.wrap_to_center(x, uc, center)
        xo = Packmol.triclinic_to_octahedral(x, uc, center)
        # the octahedral (Wigner-Seitz) image is never farther from the
        # center than the plain parallelepiped-wrapped one
        @test norm(xo - center) <= norm(xt - center) + 1e-9
        # converting the octahedral image back to the triclinic
        # representation must land exactly on the ordinary wrap of `x`
        @test Packmol.octahedral_to_triclinic(xo, uc, center) ≈ xt
    end

    # vector and PackmolSystem-argument overloads agree with the scalar form
    xs = [SVector(90.0, 5.0, 5.0), SVector(1.0, 1.0, 1.0)]
    @test Packmol.triclinic_to_octahedral(xs, uc, center) ==
          [Packmol.triclinic_to_octahedral(x, uc, center) for x in xs]

    st = structure_type(
        Packmol.src_dir * "/../test/structure_files/water.pdb";
        number=2, constraints=[InsideBox([-10.0, -10.0, -10.0], [10.0, 10.0, 10.0])],
    )
    sys = PackmolSystem([st]; output="octahedron_test.pdb", tolerance=2.0, unitcell=uc, unitcell_center=center)
    x = SVector(90.0, 5.0, 5.0)
    @test Packmol.triclinic_to_octahedral(x, sys) == Packmol.triclinic_to_octahedral(x, uc, center)
    @test Packmol.octahedral_to_triclinic(x, sys) == Packmol.octahedral_to_triclinic(x, uc, center)
end

@testitem "triclinic_to_octahedral(::PackmolSystem) / octahedral_to_triclinic(::PackmolSystem)" begin
    st = structure_type(
        Packmol.src_dir * "/../test/structure_files/water.pdb";
        number=2, constraints=[InsideBox([-10.0, -10.0, -10.0], [10.0, 10.0, 10.0])],
    )

    # No PBC: there is nothing to select a wrapping shape for.
    sys_nopbc = PackmolSystem([st]; output="octahedron_style_test.pdb", tolerance=2.0)
    @test sys_nopbc.periodic_boundary_style == :triclinic
    @test_throws ArgumentError Packmol.triclinic_to_octahedral(sys_nopbc)

    uc = Packmol.octahedral_unitcell(Float64, 100.0)
    sys = PackmolSystem([st]; output="octahedron_style_test.pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3),
    )
    @test sys.periodic_boundary_style == :triclinic

    result = Packmol.triclinic_to_octahedral(sys)
    @test result === sys
    @test sys.periodic_boundary_style == :octahedral

    result = Packmol.octahedral_to_triclinic(sys)
    @test result === sys
    @test sys.periodic_boundary_style == :triclinic
end
