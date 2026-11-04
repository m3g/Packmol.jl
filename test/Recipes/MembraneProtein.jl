@testitem "MembraneProtein setup" begin
    using Packmol
    using Unitful
    protein = joinpath(Packmol.RecipesDirectory, "test", "data", "integrin_tm.pdb")

    s = MembraneProtein(; protein, lipids = ["POPC", "POPS"], lipid_molar_ratio = [4, 1])
    @test s.protein_charge == 3 # from the residues
    @test s.margin == 15.0u"Å"
    @test s.ionic_concentration == 0.15u"mol/L"
    @test basename(s.cation_pdbfile) == "sodium.pdb"
    @test s.membrane.lipid_names == ["POPC", "POPS"]

    setup = redirect_stdout(devnull) do
        write_packmol_input(s; input = tempname() * ".inp", debug = true)
    end
    (; lipid_placements, solvent_placements, unitcell, lo, hi) = setup
    # the cell extends the margin beyond the protein in x and y
    atoms = Packmol.read_pdb(protein)
    @test lo[1] ≈ minimum(a.x for a in atoms) - 15.0
    @test hi[2] ≈ maximum(a.y for a in atoms) + 15.0
    @test unitcell[1, 1] ≈ hi[1] - lo[1]
    # lipids in both leaflets, fewer than in the cell without the protein
    upper = sum(lp.number for lp in lipid_placements if lp.upper)
    lower = sum(lp.number for lp in lipid_placements if !lp.upper)
    full = round(Int, unitcell[1, 1] * unitcell[2, 2] / ustrip(s.membrane.area_per_lipid))
    @test 0 < upper < full && 0 < lower < full
    # a neutral system
    count_of(file) = sum((sp.number for sp in solvent_placements if sp.pdbfile == file); init = 0)
    @test s.protein_charge + count_of(s.cation_pdbfile) - count_of(s.anion_pdbfile) == 0
    @test count_of(s.water_pdbfile) > 1000

    # the input file is valid for the native engine
    input = tempname() * ".inp"
    redirect_stdout(devnull) do
        write_packmol_input(s; input)
    end
    psys = Packmol.read_packmol_input(input)
    @test psys.structure_types[1].fixed.fixed
    @test psys.nmols == 1 + upper + lower + sum(sp.number for sp in solvent_placements)

    @test_throws ArgumentError MembraneProtein(; protein = "nonexistent.pdb", lipids = "POPC")
    @test_throws ArgumentError MembraneProtein(; protein, lipids = "POPC", cation_pdbfile = protein)
end

@testitem "MembraneProtein packmol" begin
    using Packmol
    protein = joinpath(Packmol.RecipesDirectory, "test", "data", "integrin_tm.pdb")
    s = MembraneProtein(; protein, lipids = "POPC", margin = 6.0, solvent_margin = 6.0)
    out = tempname() * ".pdb"
    sys = packmol(s; output = out, seed = 1, iprint = 0)
    @test sys.status == :packing_ready
    @test isfile(out)
end
