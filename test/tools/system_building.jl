@testitem "cross_section and builtin_lipid" begin
    using PDBToolsBase: read_pdb, Atom
    # One atom: a disc of radius 2 (on the grid of the computation)
    @test isapprox(Packmol.cross_section([[0.0, 0.0, 0.0]]; radius = 2.0, spacing = 0.05), 4π; rtol = 0.02)
    # Two distant atoms: two discs; outside the slab: nothing
    xs = [[0.0, 0.0, 0.0], [10.0, 0.0, 5.0]]
    @test isapprox(Packmol.cross_section(xs; spacing = 0.05), 8π; rtol = 0.02)
    @test isapprox(Packmol.cross_section(xs; z = (-1, 1), spacing = 0.05), 4π; rtol = 0.02)
    @test Packmol.cross_section(xs; z = (20, 30)) == 0.0
    # Atoms
    protein = read_pdb(joinpath(Packmol.RecipesDirectory, "test", "data", "integrin_tm.pdb"))
    @test 200 < Packmol.cross_section(protein; z = (0, 17)) < 1000
    # The built-in lipids
    popc = Packmol.builtin_lipid("popc")
    @test isfile(popc.file)
    @test (popc.head, popc.chain_ends, popc.area) == (1, [88, 131], 64.0)
    @test 25 < popc.length < 27
    for name in ("POPE", "POPS", "POPG", "DPPC", "DOPC", "DMPC")
        @test isfile(Packmol.builtin_lipid(name).file)
    end
    @test_throws ArgumentError Packmol.builtin_lipid("XYZ")
    @test Base.ispublic(Packmol, :cross_section) && Base.ispublic(Packmol, :builtin_lipid)
end
