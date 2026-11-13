@testitem "shape_descriptors and select_conformations" setup=[TestFiles] begin
    using StaticArrays: SVector
    # A rod: atoms on a line (anisotropy 1); a sphere-like set: the vertices of a cube
    rod = [SVector(0.0, 0.0, z) for z in 0.0:1.0:10.0]
    d = shape_descriptors(rod)
    @test d.anisotropy ≈ 1.0
    @test d.asphericity ≈ d.principal_variances[1]
    @test d.acylindricity ≈ 0.0 atol = 1e-12
    cube = [SVector(x, y, z) for x in (-1.0, 1.0), y in (-1.0, 1.0), z in (-1.0, 1.0)][:]
    d = shape_descriptors(cube)
    @test d.anisotropy ≈ 0.0 atol = 1e-12
    @test d.radius_of_gyration ≈ sqrt(3.0)
    disc = [SVector(cos(θ), sin(θ), 0.0) for θ in range(0, 2π; length=13)[1:12]]
    @test shape_descriptors(disc).anisotropy ≈ 0.25

    # A file with three models of water, stretched along z by different amounts
    water = TestFiles.test_file("structure_files/water.pdb")
    lines = filter(l -> startswith(l, "HETATM") || startswith(l, "ATOM"), readlines(water))
    dir = mktempdir()
    models = joinpath(dir, "models.pdb")
    open(models, "w") do io
        for (k, s) in enumerate((1.0, 3.0, 2.0))
            println(io, "MODEL     ", lpad(k, 4))
            for l in lines
                z = parse(Float64, l[47:54])
                y = parse(Float64, l[39:46])
                println(io, l[1:38], lpad(string(round(s * y; digits=3)), 8), lpad(string(round(s * z; digits=3)), 8), l[55:end])
            end
            println(io, "ENDMDL")
        end
    end
    ds = shape_descriptors(models)
    @test length(ds) == 3
    @test shape_descriptors(water) isa NamedTuple
    out = joinpath(dir, "rods.pdb")
    indices, descriptors = select_conformations(models, out; n=2)
    @test indices == [2, 3]
    @test [d.asphericity for d in descriptors] == [ds[2].asphericity, ds[3].asphericity]
    sts = structure_type(out; number=4)
    @test length(sts) == 2
    @test shape_descriptors(out)[1].asphericity ≈ ds[2].asphericity rtol = 1e-3
    indices, _ = select_conformations(models, joinpath(dir, "rod.pdb")) # n = 1 by default
    @test indices == [2]
    @test structure_type(joinpath(dir, "rod.pdb"); number=4) isa Packmol.StructureType
    indices, _ = select_conformations(models, joinpath(dir, "spheres.pdb"); n=1, shape=:sphere)
    @test length(indices) == 1
    @test_throws ArgumentError select_conformations(models, out; n=1, shape=:disc)
end
