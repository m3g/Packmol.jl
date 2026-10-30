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
        Concentration units: x
        Concentration range: 0.0 - 1.0
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
    # :cubic is the default
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 13527
    @test r1[2] == 13527
    @test r1[3] ≈ [118.81, 118.81, 118.81]u"Å"

    # The same concentration in mol/L, given with concentration_units, or with the units
    # of the system set by convert_density_table! (the density is interpolated linearly
    # in the given units, so the numbers of molecules differ slightly)
    cmol = ustrip(convert_concentration(system, 0.5, "x" => "mol/L"))
    rm(tmp_input_file, force=true)
    r2 = write_packmol_input(system; concentration = cmol, concentration_units = "mol/L", margin = 20.0, input = tmp_input_file, debug = true)
    @test all(isapprox.(r2[1:2], r1[1:2]; rtol=2e-3))
    convert_density_table!(system, "mol/L")
    rm(tmp_input_file, force=true)
    r2 = write_packmol_input(system; concentration = cmol, margin = 20.0, input = tmp_input_file, debug = true)
    @test all(isapprox.(r2[1:2], r1[1:2]; rtol=2e-3))
    convert_density_table!(system, "x")

    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true, pbc = :orthorhombic)
    @test isfile(tmp_input_file)
    @test r1[1] == 10075
    @test r1[2] == 10075
    @test r1[3] ≈ [117.37, 89.79, 118.81]u"Å"

    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true, pbc = :dodecahedral)
    @test isfile(tmp_input_file)
    @test r1[3] ≈ [118.81, 118.81, 118.81]u"Å"
    @test r1[1] + r1[2] < 13527 + 13527

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

@testitem "SolutionBoxUSC: convert_concentration" begin
    using Unitful
    test_dir = Packmol.RecipesDirectory * "/test"
    sys(table) = redirect_stderr(devnull) do
        SolutionBoxUSC(
            solute_pdbfile = "$test_dir/data/poly_h.pdb",
            solvent_pdbfile = "$test_dir/data/water.pdb",
            cossolvent_pdbfile = "$test_dir/data/water.pdb",
            density_table = table,
        )
    end
    c(x) = ustrip(x)
    # Water as solvent and cossolvent, with constant density
    system = sys(hcat(0:0.1:1, ones(11)))
    mw = 55.50824
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ x atol = 1e-3
        @test c(convert_concentration(system, x, "x" => "mol/L")) ≈ x * mw atol = 1e-2
        @test convert_concentration(system, x, "x" => "mm") ≈ x atol = 1e-3
        @test convert_concentration(system, x, "vv" => "x") ≈ x atol = 1e-3
        @test c(convert_concentration(system, x, "vv" => "mol/L")) ≈ x * mw atol = 1e-2
        @test convert_concentration(system, x, "vv" => "mm") ≈ x atol = 1e-3
        @test convert_concentration(system, x * mw, "mol/L" => "x") ≈ x atol = 1e-3
        @test convert_concentration(system, x * mw, "mol/L" => "vv") ≈ x atol = 1e-3
        @test convert_concentration(system, x * mw, "mol/L" => "mm") ≈ x atol = 1e-3
        @test convert_concentration(system, x * mw * u"mol/L", "mol/L" => "x") ≈ x atol = 1e-3
    end
    # Density linear in the molar fraction, with the same molar mass for both components
    system = sys(hcat([0.05 * i for i in 0:20], [1.0 + 0.05 * i for i in 0:20]))
    ρc, ρw = 2.0, 1.0
    M = c(system.solvent_molar_mass)
    ρ(x) = ρc * x + ρw * (1 - x)
    # Volume fraction (IUPAC): volumes of the pure components before mixing
    vv(x) = (x / ρc) / (x / ρc + (1 - x) / ρw)
    mx(x) = x / (M / (1000 * ρ(x)))
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ vv(x) atol = 1e-6
        @test c(convert_concentration(system, x, "x" => "mol/L")) ≈ mx(x) rtol = 1e-6
        @test convert_concentration(system, x, "x" => "mm") ≈ x atol = 1e-6
        @test convert_concentration(system, vv(x), "vv" => "x") ≈ x atol = 1e-6
        @test c(convert_concentration(system, vv(x), "vv" => "mol/L")) ≈ mx(x) rtol = 1e-6
        @test convert_concentration(system, vv(x), "vv" => "mm") ≈ x atol = 1e-6
        @test convert_concentration(system, mx(x), "mol/L" => "x") ≈ x atol = 1e-6
        @test convert_concentration(system, mx(x), "mol/L" => "vv") ≈ vv(x) atol = 1e-6
        @test convert_concentration(system, mx(x), "mol/L" => "mm") ≈ x atol = 1e-6
    end
    @test_throws ArgumentError convert_concentration(system, 0.5, "x" => "mol/kg")
    # convert_density_table!: the units of the system
    @test system.concentration_units == "x"
    @test get_density_table(system)[:, 1] ≈ [0.05 * i for i in 0:20]
    @test convert_density_table!(system, "mol/L") === system
    @test system.concentration_units == "mol/L"
    @test get_density_table(system)[:, 1] ≈ mx.([0.05 * i for i in 0:20]) rtol = 1e-6
    @test get_density_table(system)[:, 2] ≈ [1.0 + 0.05 * i for i in 0:20]
    convert_density_table!(system, "vv")
    @test get_density_table(system)[:, 1] ≈ vv.([0.05 * i for i in 0:20]) atol = 1e-6
    @test_throws ArgumentError convert_density_table!(system, "mol/kg")
end
