@testitem "SolutionBoxUWI" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.RecipesDirectory*"/test"

    # Default ions (SOD/CLA), neutral solute
    system = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 0)
    @test system.cation_charge == 1
    @test system.anion_charge == -1
    @test system.cation_molar_mass ≈ 22.99u"g/mol" rtol = 0.01
    @test system.anion_molar_mass ≈ 35.45u"g/mol" rtol = 0.01
    @test system.water_molar_mass ≈ 18.02u"g/mol" rtol = 0.01

    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUWI properties (Solute + Water + Ions):
    ==================================================================
        Solute pdb file: poly_h.pdb
        Solute charge: 0
        Water pdb file: HOH.pdb
        Cation pdb file and charge: SOD.pdb +1
        Anion pdb file and charge: CLA.pdb -1
        Molar masses:
            solute: 5612.80078125 g mol^-1
            water: 18.01534080505371 g mol^-1
            cation: 22.989770889282227 g mol^-1
            anion: 35.452999114990234 g mol^-1
        Ionic concentration range: 0.09953429399586157 mol L^-1 - 5.30461986641539 mol L^-1
        Density range: 1.00116 g mL^-1 - 1.19412 g mL^-1
    ==================================================================
    """

    # Invalid cation/anion charges
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", cation_charge = -1,
    )
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", anion_charge = 1,
    )
    # Custom ion PDB file without an explicit charge
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", cation_pdbfile = "$test_dir/data/SOD.pdb",
    )

    # Neutral solute, no added salt: only water fills the box
    system = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 0)
    tmp_input_file = tempname()*".inp"
    rm(tmp_input_file, force=true)
    r = write_packmol_input(system; ionic_concentration = 0.0, margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r[1] > 0  # water
    @test r[2] == 0  # cations
    @test r[3] == 0  # anions
    @test r[4] ≈ [118.81, 118.81, 118.81]u"Å"
    input_text = read(tmp_input_file, String)
    @test occursin("structure $(Packmol.RecipesDirectory)/test/data/poly_h.pdb", input_text)
    @test !occursin("structure $(system.cation_pdbfile)", input_text)
    @test !occursin("structure $(system.anion_pdbfile)", input_text)

    # Neutral solute, physiological ionic strength: equal numbers of cations and anions
    rm(tmp_input_file, force=true)
    r = write_packmol_input(system; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r[2] == r[3]
    @test r[2] > 0
    input_text = read(tmp_input_file, String)
    @test occursin("structure $(system.cation_pdbfile)", input_text)
    @test occursin("structure $(system.anion_pdbfile)", input_text)

    # The generated input file must be valid Packmol syntax for the native engine
    # itself, not just the (more lenient) legacy Fortran binary — regression test
    # for a stray comma once present in the `pbc` line, which the native parser
    # rejected outright.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == 1 + r[1] + r[2] + r[3]
    rm(tmp_input_file, force=true)

    # Charged solute (+5): extra anions neutralize it, on top of the bulk salt
    system_charged = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 5)
    rm(tmp_input_file, force=true)
    r_charged = write_packmol_input(system_charged; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test r_charged[3] - r_charged[2] == 5  # 5 extra anions relative to cations

    # Negatively charged solute (-3): extra cations neutralize it
    system_neg = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = -3)
    rm(tmp_input_file, force=true)
    r_neg = write_packmol_input(system_neg; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test r_neg[2] - r_neg[3] == 3  # 3 extra cations relative to anions

    # Custom divalent cation (charge +2) neutralizing a negative solute charge: an odd
    # magnitude cannot be exactly cancelled by ions of charge magnitude 2
    system_divalent = SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solute_charge = -3,
        cation_pdbfile = "$test_dir/data/SOD.pdb",
        cation_charge = 2,
    )
    @test_throws ArgumentError write_packmol_input(
        system_divalent; ionic_concentration = 0.15u"mol/L", margin = 20.0, input = tmp_input_file,
    )

    rm(tmp_input_file, force=true)
end
