@testitem "Consistency tests" begin
    using ShowMethodTesting
    using PDBTools
    using Packmol

    test_dir = Packmol.PackmolInputCreatorDirectory*"/test"

    # system with water only, with constant density
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat(0:0.1:1, ones(11))
    )
    mw = 55.508250191225926
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ x atol = 1e-3
        @test convert_concentration(system, x, "x" => "mol/L") ≈ x * mw atol = 1e-3
        @test convert_concentration(system, x, "x" => "mm") ≈ x atol = 1e-3

        @test convert_concentration(system, x, "vv" => "x") ≈ x atol = 1e-3
        @test convert_concentration(system, x, "vv" => "mol/L") ≈ x * mw atol = 1e-3
        @test convert_concentration(system, x, "vv" => "mm") ≈ x atol = 1e-3

        @test convert_concentration(system, x * mw, "mol/L" => "x") ≈ x atol = 1e-3 
        @test convert_concentration(system, x * mw, "mol/L" => "vv") ≈ x atol = 1e-3
        @test convert_concentration(system, x * mw, "mol/L" => "mm") ≈ x atol = 1e-3
    end

    # system with ideal solution 
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat([0.0 + 0.05*i for i in 0:20], [1.0 + 0.05*i for i in 0:20])
    )
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.cossolvent_molar_mass ≈ 18.01 atol = 0.01
    @test density_pure_solvent(system) ≈ 1.0
    @test density_pure_cossolvent(system) ≈ 2.0

    # Concentration conversions in this ideal system, from the molar fraction, x
    ρc = density_pure_cossolvent(system) # g / mL
    ρw = density_pure_solvent(system) # g / mL
    M = system.solvent_molar_mass # g / mol
    ρ(x) = ρc*x + ρw*(1-x) # g / mL
    vv(x) = (x / ρc) / (x / ρc + (1-x) / ρw) 
    v(x) = M / (1000*ρ(x)) # L/mol: volume of 1 mol (c + w) of solution
    mx(x) = x / v(x) # molarity of cossolute
    mm(x) = x # molality of cossolute

    # Concentration conversions
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ vv(x) 
        @test convert_concentration(system, x, "x" => "mol/L") ≈ mx(x)
        @test convert_concentration(system, x, "x" => "mm") ≈ x

        @test convert_concentration(system, vv(x), "vv" => "x") ≈ x 
        @test convert_concentration(system, vv(x), "vv" => "mol/L") ≈ mx(x) 
        @test convert_concentration(system, vv(x), "vv" => "mm") ≈ x

        @test convert_concentration(system, mx(x), "mol/L" => "x") ≈ x 
        @test convert_concentration(system, mx(x), "mol/L" => "vv") ≈ vv(x) 
        @test convert_concentration(system, mx(x), "mol/L" => "mm") ≈ x
    end

    # Water as cossolvent
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
    system = @test_logs (:warn,) SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw
    )
    @test convert_concentration(system, 1.0, "x" => "mol/L") ≈ 55.40278451586260
    dw[end, 1] = 0.99
    @test_throws ArgumentError SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw
    )
    dw[end, 1] = 1.0
    dw[begin, 1] = 0.1
    @test_throws ArgumentError SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw
    )
    dw[begin, 1] = 0.0
    convert_density_table!(system, "mol/L")
    dw = copy(system.density_table)
    @test_throws ArgumentError SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw,
        concentration_units = "x",
    )
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw,
        concentration_units = "mol/L",
    )
    @test density_pure_cossolvent(system) ≈ 0.9981
    dw[end, 1] = dw[end, 1] + 0.01
    @test_throws ArgumentError SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/ethanol.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = dw,
        concentration_units = "mol/L",
    )

    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat(0:0.1:1, ones(11)),
        concentration_units = "x",
    )
    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUSC properties (Solute + Solvent + Cossolvent):
    ==================================================================
    Solute pdb file: poly_h.pdb
    Solvent pdb file: water.pdb
    Cossolvent pdb file: water.pdb
    Density of pure solvent: 1.0 g/mL
    Density of pure cossolvent: 1.0 g/mL
    Molar masses: 
        solute: 5612.791939999981 g/mol
        solvent: 18.01534 g/mol
        cossolvent: 18.01534 g/mol
    Concentration units: x (molar fraction)
    Cocentration range: 0.0 - 1.0
    ==================================================================
    """

end

@testitem "Write packmol input" begin
    using PDBTools
    using MolSimToolkit.PackmolInputCreator
    using DelimitedFiles
    test_dir = PackmolInputCreator.PackmolInputCreatorDirectory*"/test"

    # Ethanol-water mixture
    density_table = readdlm("$test_dir/data/water_ethanol.dat", comments=true, comment_char='#')
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
        density_table = copy(density_table),
    )
    Mc = 1000 * density_pure_cossolvent(system) / system.cossolvent_molar_mass # mol / L pure ethanol

    # Test concentration conversions for real data
    @test convert_concentration(system, 1.0, "x" => "mol/L") ≈ Mc

    mm = 0.3997224931406948
    vv = 0.45671854335897716
    x = 0.2066
    M = 8.12907194485856
    ρ = 0.9369

    @test convert_concentration(system, x, "x" => "vv") ≈ vv
    @test convert_concentration(system, x, "x" => "mol/L") ≈ M
    @test convert_concentration(system, x, "x" => "mm") ≈ mm 

    @test convert_concentration(system, vv, "vv" => "x") ≈ x
    @test convert_concentration(system, vv, "vv" => "mol/L") ≈ M
    @test convert_concentration(system, vv, "vv" => "mm") ≈ mm 

    @test convert_concentration(system, M, "mol/L" => "x") ≈ x
    @test convert_concentration(system, M, "mol/L" => "mm") ≈ mm
    @test convert_concentration(system, M, "mol/L" => "vv") ≈ vv

    @test convert_concentration(system, mm, "mm" => "x") ≈ x
    @test convert_concentration(system, mm, "mm" => "mol/L") ≈ M
    @test convert_concentration(system, mm, "mm" => "vv") ≈ vv

    tmp_input_file = tempname()

    convert_density_table!(system, "x")
    r1 = write_packmol_input(system; concentration = 0.0, box_sides=[120,120,120], input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 57322
    
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true, cubic = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 13531
    @test r1[2] == 13531
    @test r1[3] ≈ [118.81, 118.81, 118.81]

    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r1[1] == 10080
    @test r1[2] == 10080
    @test r1[3] ≈ [117.37, 89.79, 118.81]

    convert_density_table!(system, "mol/L")
    r2 = write_packmol_input(system; concentration = 13.488667939471432, margin = 20.0, input = tmp_input_file, debug = true)
    @test all(isapprox.(r2,r1,rtol=0.005))

    convert_density_table!(system, "vv")
    r3 = write_packmol_input(system; concentration = 0.7635032204047275, margin = 20.0, input = tmp_input_file, debug = true)
    @test all(isapprox.(r3,r1,rtol=0.005))

    convert_density_table!(system, "mm")
    r4 = write_packmol_input(system; concentration = 0.7188817400010237, margin = 20.0, input = tmp_input_file, debug = true)
    @test all(isapprox.(r4,r1,rtol=0.005))

    convert_density_table!(system, "mol/L")
    system.concentration_units = "x"
    @test_throws ArgumentError write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true)
    system.concentration_units = "mol/L"
    system.density_table .= density_table
    @test_throws ArgumentError write_packmol_input(system; concentration = 0.5, margin = 20.0, input = tmp_input_file, debug = true)


end