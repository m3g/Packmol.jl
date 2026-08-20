@testitem "SolutionBoxUSC" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.RecipesDirectory*"/test"

    mw = 55.508250191225926
    # system with ideal solution 
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat([0.0 + 0.05*i for i in 0:20], [1.0 + 0.05*i for i in 0:20])
    )
    @test system.solvent_molar_mass ≈ 18.01u"g/mol" rtol=0.01
    @test system.cossolvent_molar_mass ≈ 18.01u"g/mol" rtol = 0.01
    @test density_pure_solvent(system) ≈ 1.0u"g/mL"


    # System of water in water, easy to test
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat(0:0.1:1, ones(11)),
        concentration_units = "x",
    )

    @test Packmol.density_pure_solvent(system) ≈ 1.0u"g/mL" atol=0.01u"g/mL"
    @test Packmol.density_pure_cossolvent(system) ≈ 1.0u"g/mL" atol=0.01u"g/mL"

    # missing pure solute and cossolvent densities
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat(0.1:0.1:0.9, ones(9)),
        concentration_units = "x",
    )
    @test ismissing(Packmol.density_pure_solvent(system))
    @test ismissing(Packmol.density_pure_cossolvent(system))

    # Ethanol/Water, Water as cossolvent
    dw = [
        0.0     0.7906
        0.2214  0.8195
        0.3902  0.845
        0.5231  0.8685
        0.6305  0.8923
        0.7191  0.9151
        0.7934  0.9369
        0.8566  0.9537
        0.911   0.9685
        0.9584  0.982
        1.0     0.9981
    ]
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw
    )
    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUSC properties (Solute + Solvent + Cossolvent):
    ==================================================================
        Solute pdb file: poly_h.pdb
        Solvent pdb file: ethanol.pdb
        Cossolvent pdb file: water.pdb
        Molar masses:
            solute: 5612.801f0 g mol^-1
            solvent: 46.069218f0 g mol^-1
            cossolvent: 18.01534f0 g mol^-1
        Concentration range: 0.0 mol L^-1 - 55.40278204007167 mol L^-1
        Density range: 0.7906 g mL^-1 - 0.9981 g mL^-1
    ==================================================================
    """

    @test first(system.density_table.concentration) ≈ 0.0u"mol/L" atol=0.01u"mol/L"
    @test last(system.density_table.concentration) ≈ 55.40u"mol/L" rtol=0.01
    @test parse_show(system.density_table) ≈ """
        ==================================================================
        Density table:
        ==================================================================
        Concentration (mol L^-1) |      Density (g mL^-1)
                  0.000          |           0.791
                  4.552          |           0.820
                  9.388          |           0.845
                 14.471          |           0.869
                 19.823          |           0.892
                 25.412          |           0.915
                 31.218          |           0.937
                 37.069          |           0.954
                 43.014          |           0.969
                 49.063          |           0.982
                 55.403          |           0.998
        ==================================================================
    """
    @test Packmol.interpolate_density(system, 0.0u"mol/L", "mol/L") ≈ 0.791u"g/mL" rtol=0.01
    @test Packmol.interpolate_density(system, 55.402u"mol/L", "mol/L") ≈ 0.9981u"g/mL" rtol=0.01
    @test Packmol.interpolate_density(system, 25.412u"mol/L", "mol/L") ≈ 0.915u"g/mL" rtol=0.01
    @test Packmol.interpolate_density(system, 0.0, "x") ≈ 0.791u"g/mL" atol=0.01u"g/mL"
    @test Packmol.interpolate_density(system, 1.0, "x") ≈ 0.9981u"g/mL" atol=0.01u"g/mL"
    @test Packmol.interpolate_density(system, 0.5, "x") ≈ 0.863u"g/mL" atol=0.01u"g/mL"

    tmp_input_file = tempname() * ".inp"
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true, cubic = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 13527
    @test r1[2] == 13527
    @test r1[3] ≈ [118.81, 118.81, 118.81]u"Å"

    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 10075
    @test r1[2] == 10075
    @test r1[3] ≈ [117.37, 89.79, 118.81]u"Å"

    # The generated input file must be valid Packmol syntax for the native engine
    # itself, not just the (more lenient) legacy Fortran binary — regression test
    # for a stray comma once present in the `pbc` line, which the native parser
    # rejected outright.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == 1 + r1[1] + r1[2]
    rm(tmp_input_file, force=true)

    # In this test, if we (incorrectly) provide the concentration in mol/L,
    # the conversion of the density table will not change the values
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw,
        concentration_units = "mol/L",
    )
    @test first(system.density_table.concentration) ≈ 0.0u"mol/L" atol=0.01u"mol/L"
    @test last(system.density_table.concentration) ≈ 1.0u"mol/L" rtol=0.01

end
