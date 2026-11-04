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
