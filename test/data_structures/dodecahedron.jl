@testitem "dodecahedral_unitcell" begin
    uc = Packmol.dodecahedral_unitcell(Float64, 100.0)
    a, b, c, α, β, γ = Packmol._unitcell_abc_angles(uc)
    @test a ≈ 100.0
    @test b ≈ 100.0
    @test c ≈ 100.0
    @test α ≈ 60.0
    @test β ≈ 60.0
    @test γ ≈ 90.0
end

@testitem "parse_pbc_dodecahedral" begin
    using StaticArrays: SVector
    uc, center = Packmol.parse_pbc_dodecahedral(Float64, ["1.0", "2.0", "3.0", "100.0"], 3)
    @test uc ≈ Packmol.dodecahedral_unitcell(Float64, 100.0)
    @test center ≈ SVector(1.0, 2.0, 3.0)
    @test_throws ArgumentError Packmol.parse_pbc_dodecahedral(Float64, ["1.0", "2.0", "3.0", "100.0"], 2)
    @test_throws ArgumentError Packmol.parse_pbc_dodecahedral(Float64, ["1.0", "2.0", "3.0"], 3)
end

@testitem "pbc dodecahedral input keyword" begin
    using Packmol: read_packmol_input
    dir = dirname(Packmol.src_dir * "/../test/input_files/water_box.inp")
    original = read(Packmol.src_dir * "/../test/input_files/water_box.inp", String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))
        file = joinpath(tmp, "dodeca.inp")
        write(file, original * "\npbc dodecahedral 0. 0. 0. 100.\n")
        sys = read_packmol_input(file)
        @test !isnothing(sys.unitcell)
        a, b, c, α, β, γ = Packmol._unitcell_abc_angles(sys.unitcell)
        @test a ≈ b ≈ c ≈ 100.0
        @test α ≈ β ≈ 60.0
        @test γ ≈ 90.0
        @test sys.unitcell_center == zeros(3)
    end
end

@testitem "triclinic_to_dodecahedral / dodecahedral_to_triclinic round trip" begin
    using StaticArrays: SVector
    using LinearAlgebra: norm
    using Random

    d = 100.0
    uc = Packmol.dodecahedral_unitcell(Float64, d)
    center = SVector(0.0, 0.0, 0.0)

    Random.seed!(1)
    for _ in 1:2000
        x = SVector(d * (rand() - 0.5), d * (rand() - 0.5), d * (rand() - 0.5))
        xt = Packmol.wrap_to_center(x, uc, center)
        xd = Packmol.triclinic_to_dodecahedral(x, uc, center)
        # the dodecahedral (Wigner-Seitz) image is never farther from the
        # center than the plain parallelepiped-wrapped one
        @test norm(xd - center) <= norm(xt - center) + 1e-9
        # converting the dodecahedral image back to the triclinic
        # representation must land exactly on the ordinary wrap of `x`
        @test Packmol.dodecahedral_to_triclinic(xd, uc, center) ≈ xt
    end

    # vector and PackmolSystem-argument overloads agree with the scalar form
    xs = [SVector(90.0, 5.0, 5.0), SVector(1.0, 1.0, 1.0)]
    @test Packmol.triclinic_to_dodecahedral(xs, uc, center) ==
          [Packmol.triclinic_to_dodecahedral(x, uc, center) for x in xs]

    st = structure_type(
        Packmol.src_dir * "/../test/structure_files/water.pdb";
        number=2, constraints=[InsideBox([-10.0, -10.0, -10.0], [10.0, 10.0, 10.0])],
    )
    sys = PackmolSystem([st]; output="dodecahedron_test.pdb", tolerance=2.0, unitcell=uc, unitcell_center=center)
    x = SVector(90.0, 5.0, 5.0)
    @test Packmol.triclinic_to_dodecahedral(x, sys) == Packmol.triclinic_to_dodecahedral(x, uc, center)
    @test Packmol.dodecahedral_to_triclinic(x, sys) == Packmol.dodecahedral_to_triclinic(x, uc, center)
end
