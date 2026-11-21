@testitem "scaled constraints" begin
    using Packmol
    using Packmol: _scaled_constraint, constraint_penalty, InsideSphere, OutsideSphere, InsideBox,
        InsideCube, InsideCylinder, InsideEllipsoid, AbovePlane, BelowPlane
    using StaticArrays
    import Random
    rng = Random.Xoshiro(1)
    o = SVector(1.0, -2.0, 0.5)
    s = 0.8
    constraints = (
        InsideSphere([1.0, 2.0, 3.0], 5.0), OutsideSphere([0.0, 0.0, 0.0], 4.0),
        InsideBox([-3.0, -2.0, -1.0], [4.0, 5.0, 6.0]), InsideCube([-2.0, -2.0, -2.0], 5.0),
        InsideCylinder([0.0, 0.0, -3.0], [1.0, 1.0, 1.0], 3.0, 8.0),
        InsideEllipsoid([1.0, 0.0, 0.0], 3.0, 4.0, 5.0, 1.0),
        AbovePlane([1.0, 2.0, 3.0], 2.0), BelowPlane([0.0, 0.0, 1.0], -1.0),
    )
    for c in constraints
        cs = _scaled_constraint(c, s, o)
        # x satisfies c if and only if its image satisfies the scaled constraint
        for _ in 1:2000
            x = SVector{3}(10 .* (2 .* rand(rng, 3) .- 1))
            @test (constraint_penalty(c, x) == 0) == (constraint_penalty(cs, o + s * (x - o)) == 0)
        end
    end
    @test_throws ArgumentError _scaled_constraint(Packmol.AboveGaussian([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 0.0, 1.0, 0.0, 1.0), s, o)
end

@testitem "packing by compression" setup=[TestFiles] begin
    using Packmol
    using Packmol: read_packmol_input, InsideSphere
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    cp(water, joinpath(dir, "water.pdb"))
    input = joinpath(dir, "use_compression.inp")
    write(input, """
        tolerance 2.0
        output compress_out.pdb
        use_compression 0.5
        structure water.pdb
          number 300
          inside sphere 0. 0. 0. 14.
        end structure
        """)
    sys = read_packmol_input(input)
    @test sys.use_compression == (0.5, 0.05)
    sys = redirect_stdout(devnull) do
        packmol(sys)
    end
    @test sys.status == :packing_ready
    # the geometry of the input
    @test sys.structure_types[1].constraints[1].radius ≈ 14.0
    @test isfile(joinpath(dir, "compress_out.pdb"))
    # periodic cell, and from the Julia API
    sys = redirect_stdout(devnull) do
        packmol(input; use_compression = (0.7, 0.1), pbc = [-15, -15, -15, 15, 15, 15])
    end
    @test sys.status == :packing_ready
    @test sys.unitcell ≈ [30.0 0 0; 0 30.0 0; 0 0 30.0]
    # the defaults: `use_compression` alone, or `use_compression = true`
    write(input, "output x.pdb\nuse_compression\nstructure water.pdb\n number 1\nend structure\n")
    @test read_packmol_input(input).use_compression == (0.5, 0.05)
    st = structure_type(joinpath(dir, "water.pdb"); number = 1)
    @test PackmolSystem([st]; output = "x.pdb", use_compression = true).use_compression == (0.5, 0.05)
    # errors
    @test_throws ArgumentError packmol(input; use_compression = 1.5)
    write(input, "output x.pdb\nuse_compression 0.5 2.0\nstructure water.pdb\n number 1\nend structure\n")
    @test_throws ArgumentError read_packmol_input(input)
    # not with voxel_place (an error before the voxel placement)
    write(input, "output x.pdb\nuse_compression\nvoxel_place\nstructure water.pdb\n number 1\nend structure\n")
    @test_throws "not supported with voxel_place" redirect_stdout(devnull) do
        packmol(input)
    end
end

@testitem "compression with fixed molecules, dodecahedral cell" setup=[TestFiles] begin
    using Packmol
    using LinearAlgebra: norm
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    cp(water, joinpath(dir, "water.pdb"))
    input = joinpath(dir, "compress_fixed.inp")
    write(input, """
        tolerance 2.0
        output compress_fixed_out.pdb
        pbc dodecahedral 0. 0. 0. 30.
        use_compression
        structure water.pdb
          number 1
          center
          fixed 0. 0. 0. 0. 0. 0.
        end structure
        structure water.pdb
          number 500
        end structure
        """)
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    @test sys.unitcell ≈ Packmol.dodecahedral_unitcell(Float64, 30.0)
    atoms = get_atoms(sys)
    # the fixed molecule where it was placed (scaled with the geometry during the
    # compression, and restored exactly), and the others at the tolerance from it
    fixed = [[a.x, a.y, a.z] for a in atoms[1:3]]
    @test norm(sum(fixed) / 3) < 1.0
    reference = read_packmol_input(input).structure_types[1].reference_coordinates
    @test sys.structure_types[1].reference_coordinates == reference
    @test all(fixed[i] ≈ reference[i] for i in 1:3)
    @test all(sys.atoms[i].radius == 1.0 for i in 1:3)
    @test minimum(norm(x - [a.x, a.y, a.z]) for x in fixed for a in atoms[4:end]) > 1.9
end
