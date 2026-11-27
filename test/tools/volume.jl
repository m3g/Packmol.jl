@testitem "guess_volume" begin
    using Packmol
    water = joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb")
    # 1000 molecules of 18.015 g/mol at 1 g/mL: 29915 Å³
    g = redirect_stdout(devnull) do
        guess_volume([18.015 => 1000])
    end
    @test g.mass ≈ 18015.0
    @test g.volume ≈ 18015.0 / 6.02214076e23 * 1e24 rtol = 1e-10
    @test g.side ≈ cbrt(g.volume)
    # half the volume at twice the density, the molar mass read from the file
    g2 = redirect_stdout(devnull) do
        guess_volume([water => 1000]; density = 2.0)
    end
    @test g2.volume ≈ g.volume / 2 rtol = 1e-3
    # from a system
    sys = PackmolSystem([structure_type(water; number = 1000, constraints = [InsideCube([0., 0., 0.], 40.)])];
        tolerance = 2.0, output = tempname() * ".pdb")
    g3 = redirect_stdout(devnull) do
        guess_volume(sys)
    end
    @test g3.volume ≈ g2.volume * 2 rtol = 1e-10
    @test_throws ArgumentError guess_volume([18.0 => 10]; density = 0.0)
end

@testitem "free_volume" begin
    using Packmol
    using StaticArrays
    water = joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb")
    out = tempname() * ".pdb"
    quiet(f) = redirect_stdout(f, devnull)

    # A sphere of radius 10 (4188.8 Å³), and a spherical shell between 5 and 10
    inside = structure_type(water; number = 100, constraints = [InsideSphere([0., 0., 0.], 10.)])
    shell = structure_type(water; number = 10,
        constraints = [InsideSphere([0., 0., 0.], 10.), OutsideSphere([0., 0., 0.], 5.)])
    r = quiet(() -> free_volume(PackmolSystem([inside, shell]; tolerance = 2.0, output = out); npoints = 400_000))
    @test abs(r[1].volume - 4π / 3 * 10^3) < 4 * r[1].error
    @test abs(r[2].volume - 4π / 3 * (10^3 - 5^3)) < 4 * r[2].error
    @test r[1].density ≈ 100 * 18.015 / 6.02214076e23 / (r[1].volume * 1e-24) rtol = 1e-2
    @test all(isempty(x.note) for x in r)

    # A fixed atom at the center removes a sphere of the radius of the tolerance
    atom = joinpath(Packmol.RecipesDirectory, "test", "data", "SOD.pdb")
    fixed = structure_type(atom; number = 1, fixed = ([0., 0., 0.], [0., 0., 0.]))
    r = quiet(() -> free_volume(PackmolSystem([fixed, inside]; tolerance = 2.0, output = out); npoints = 400_000))
    @test length(r) == 1
    @test abs(r[1].volume - 4π / 3 * (10^3 - 2^3)) < 4 * r[1].error

    # The periodic cell, with a fixed atom on a face (its periodic image counts too)
    free = structure_type(water; number = 10)
    face = structure_type(atom; number = 1, fixed = ([0., 0., 9.5], [0., 0., 0.]))
    cell = [20. 0. 0.; 0. 20. 0.; 0. 0. 20.]
    r = quiet(() -> free_volume(PackmolSystem([face, free]; tolerance = 2.0, output = out, unitcell = cell);
        npoints = 400_000))
    @test abs(r[1].volume - (8000 - 4π / 3 * 2^3)) < 4 * r[1].error

    # No estimate: an unbounded region, and constraints of some atoms only
    plane = structure_type(water; number = 10, constraints = [AbovePlane([0., 0., 1.], 0.)])
    some = structure_type(water; number = 10, atom_constraints = [[1] => InsideSphere([0., 0., 0.], 5.)])
    r = quiet(() -> free_volume(PackmolSystem([plane, some]; tolerance = 2.0, output = out); npoints = 1000))
    @test isnan(r[1].volume) && contains(r[1].note, "not bounded")
    @test isnan(r[2].volume)
    r = quiet(() -> free_volume(PackmolSystem([some]; tolerance = 2.0, output = out, unitcell = cell); npoints = 1000))
    @test isnan(r[1].volume) && contains(r[1].note, "some atoms only")
    # some atoms, and the whole molecule: an upper bound
    both = structure_type(water; number = 10, constraints = [InsideSphere([0., 0., 0.], 10.)],
        atom_constraints = [[1] => InsideSphere([0., 0., 0.], 5.)])
    r = quiet(() -> free_volume(PackmolSystem([both]; tolerance = 2.0, output = out); npoints = 100_000))
    @test contains(r[1].note, "upper bound")
    @test abs(r[1].volume - 4π / 3 * 10^3) < 4 * r[1].error
end

@testitem "JavaScript volume tools agree with free_volume" begin
    # docs/src/assets/volume_tools.js (the online calculators of the docs) against the
    # Julia estimates, on the inputs of the tutorials (skipped without Node.js)
    node = Sys.which("node")
    if isnothing(node)
        @info "Node.js not found: the JavaScript volume tools are not tested"
    else
        docs = joinpath(pkgdir(Packmol), "docs")
        reference = tempname() * ".json"
        run(pipeline(`$(Base.julia_cmd()) --project=$(pkgdir(Packmol)) $(joinpath(docs, "test_volume_tools.jl")) $reference`; stdout=devnull))
        # JSON has no NaN (no estimate): null
        write(reference, replace(read(reference, String), "\"NaN\"" => "null"))
        out = IOBuffer()
        p = run(pipeline(ignorestatus(`$node $(joinpath(docs, "test_volume_tools.js")) $reference`); stdout=out, stderr=out))
        result = String(take!(out))
        p.exitcode == 0 || println(result)
        @test p.exitcode == 0
        @test occursin("all tests passed", result)
    end
end

@testitem "molecular_volume" begin
    using Packmol
    import PDBToolsBase
    using PDBToolsBase: Atom
    quiet(f) = redirect_stdout(f, devnull)
    # A single atom: the volume of its sphere (a carbon, of radius 1.7 Å), with any probe
    c = [Atom(name = "C", x = 0.0f0, y = 0.0f0, z = 0.0f0)]
    @test molecular_volume(c; probe = 0.0, resolution = 0.1) ≈ 4π / 3 * 1.7^3 rtol = 0.01
    @test molecular_volume(c; resolution = 0.1) ≈ 4π / 3 * 1.7^3 rtol = 0.01
    # Two atoms in contact: with the probe, the gap between them is part of the volume
    cc = [Atom(name = "C", x = 0.0f0, y = 0.0f0, z = 0.0f0), Atom(name = "C", x = 3.4f0, y = 0.0f0, z = 0.0f0)]
    vdw = molecular_volume(cc; probe = 0.0, resolution = 0.1)
    @test vdw ≈ 2 * 4π / 3 * 1.7^3 rtol = 0.01
    @test molecular_volume(cc; resolution = 0.1) > 1.05 * vdw
    # A ring of atoms around a cavity where the probe fits: the cavity is not part of the
    # volume, with or without the probe (vdW), in the limit of the excluded gaps
    ring = [Atom(name = "C", x = Float32(5cos(θ)), y = Float32(5sin(θ)), z = 0.0f0) for θ in range(0, 2π; length = 13)[1:12]]
    @test molecular_volume(ring; resolution = 0.2) < 1.5 * molecular_volume(ring; probe = 0.0, resolution = 0.2)
    # A protein: about its volume in solution (0.73 mL/g), and converging with the resolution:
    # with the probes on the creases of the accessible surface, the default (0.4 Å) is within
    # 2% of the volume at 0.2 Å (without them, 4-5% larger than the converged one)
    protein = PDBToolsBase.read_pdb(joinpath(Packmol.RecipesDirectory, "test", "data", "poly_h.pdb"))
    v4 = molecular_volume(protein)
    v2 = molecular_volume(protein; resolution = 0.2)
    @test v4 ≈ v2 rtol = 0.02
    @test v4 > v2
    @test v2 ≈ 0.73 * PDBToolsBase.mass(protein) / 0.602214076 rtol = 0.25
    # Of a system and of structure files with numbers
    water = joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb")
    vw = molecular_volume(PDBToolsBase.read_pdb(water))
    r = quiet(() -> molecular_volume([water => 10]))
    @test r.volume ≈ 10 * vw
    @test only(r.molecules).volume ≈ vw
    sys = PackmolSystem([structure_type(water; number = 10, constraints = [InsideCube([0., 0., 0.], 20.)])];
        tolerance = 2.0, output = tempname() * ".pdb")
    @test quiet(() -> molecular_volume(sys)).volume ≈ 10 * vw
    @test_throws ArgumentError molecular_volume(c; resolution = 0.0)
    @test_throws ArgumentError molecular_volume(c; probe = -1.0)
    @test_throws ArgumentError molecular_volume([Atom(name = "X")])
end
