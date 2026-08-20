@testitem "SolutionBoxUS" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.RecipesDirectory*"/test"

    # system with water only
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=55.5u"mol/L",
    ) 
    @test system.solvent_molar_mass ≈ 18.01u"g/mol" rtol = 0.01
    @test system.density ≈ 1.0u"g/mL" rtol = 0.01
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0u"g/mL",
    ) 
    @test system.solvent_molar_mass ≈ 18.01u"g/mol" rtol = 0.01
    @test system.density ≈ 1.0u"g/mL" rtol = 0.01
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
    ) 
    @test system.solvent_molar_mass ≈ 18.01u"g/mol" rtol = 0.01
    @test system.density ≈ 1.0u"g/mL" rtol = 0.01
    @test_throws ArgumentError SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=0.0,
    ) 
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=55.5u"mol/L",
        solvent_molar_mass=18.01534,
        solute_molar_mass=5612.79194,
    ) 
    @test system.density ≈ 1.0u"g/mL" rtol = 0.01

    # Test show methods
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
    ) 
    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUS properties (Solute + Solvent):
    ==================================================================
        Solute pdb file: poly_h.pdb
        Solvent pdb file: water.pdb
        Density of pure solvent: 1.0 g mL^-1
        Molarity of pure solvent: 55.508250191225926 mol L^-1
        Molar masses: 
            solute: 5612.791939999981 g mol^-1
            solvent: 18.01534 g mol^-1
    ==================================================================
    """

    # System with water only
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
    ) 
    tmp_input_file = tempname()*".inp"
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; margin = 20.0, input = tmp_input_file, debug = true)
    @test r1[1] == 41543
    @test r1[2] ≈ [117.37, 89.79, 118.81]u"Å" 
    @test isfile(tmp_input_file)
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; margin = 20.0, input = tmp_input_file, debug = true, cubic = true)
    @test r1[1] == 55750
    @test r1[2] ≈ [118.81, 118.81, 118.81]u"Å"
    @test isfile(tmp_input_file)
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; margin = 2.0u"nm", input = tmp_input_file, debug = true, cubic = true)
    @test r1[1] == 55750
    @test r1[2] ≈ [118.81, 118.81, 118.81]u"Å"
    @test isfile(tmp_input_file)

    # The generated input file must be valid Packmol syntax for the native engine
    # itself, not just the (more lenient) legacy Fortran binary — regression test
    # for a stray comma once present in the `pbc` line, which the native parser
    # rejected outright.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == 1 + r1[1]
    rm(tmp_input_file, force=true)

    # packmol(::SolutionBoxUS) builds and packs the system directly: no input
    # file is left behind, and the packed output file is written.
    small_system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/octane.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density = 1.0,
    )
    tmp_output_file = tempname()*".pdb"
    packmol(small_system; box_sides = [25.0, 25.0, 25.0], output = tmp_output_file, iprint = 1000, nloop = 20)
    @test isfile(tmp_output_file)
    rm(tmp_output_file, force=true)

end
