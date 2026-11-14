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
    r1 = write_packmol_input(system; margin = 20.0, input = tmp_input_file, debug = true, pbc = :orthorhombic)
    @test r1[1] == 41543
    @test r1[2] ≈ [117.37, 89.79, 118.81]u"Å"
    @test isfile(tmp_input_file)
    rm(tmp_input_file, force=true)
    # :cubic is the default
    r1 = write_packmol_input(system; margin = 20.0, input = tmp_input_file, debug = true)
    @test r1[1] == 55750
    @test r1[2] ≈ [118.81, 118.81, 118.81]u"Å"
    @test isfile(tmp_input_file)
    rm(tmp_input_file, force=true)
    r1 = write_packmol_input(system; margin = 2.0u"nm", input = tmp_input_file, debug = true, pbc = :cubic)
    @test r1[1] == 55750
    @test r1[2] ≈ [118.81, 118.81, 118.81]u"Å"
    @test isfile(tmp_input_file)
    rm(tmp_input_file, force=true)
    # :dodecahedral: the smallest cell with every atom of the (elongated) solute at least
    # the margin from its faces, larger than the cube of the box of the solute (118.81 Å),
    # whose corners are outside a dodecahedron of that size
    r1 = write_packmol_input(system; margin = 20.0, input = tmp_input_file, debug = true, pbc = :dodecahedral)
    @test r1[2] ≈ [150.574, 150.574, 150.574]u"Å" rtol = 1e-4
    @test isfile(tmp_input_file)

    # The generated input file must be valid Packmol syntax for the native engine
    # itself, not just the (more lenient) legacy Fortran binary — regression test
    # for a stray comma once present in the `pbc` line, which the native parser
    # rejected outright.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == 1 + r1[1]
    a, b, c, α, β, γ = Packmol._unitcell_abc_angles(psys.unitcell)
    @test isapprox(a, 150.574; rtol = 1e-4) && b ≈ a && c ≈ a
    Packmol._check_fixed_inside_cell(psys) # the solute is inside the cell (no error)
    @test α ≈ β ≈ 60.0
    @test γ ≈ 90.0
    rm(tmp_input_file, force=true)

    # packmol(::SolutionBoxUS) builds and packs the system directly: no input
    # file is left behind, and the packed output file is written.
    small_system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/octane.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density = 1.0,
    )
    tmp_output_file = tempname()*".pdb"
    # packmol(::Recipe) forwards packmol(::PackmolSystem)'s own return value:
    # the built PackmolSystem itself (with its packing outcome in `.status`),
    # not a bare Bool, so the caller can inspect/reuse the packed system.
    packed_sys = packmol(small_system; box_sides = [25.0, 25.0, 25.0], output = tmp_output_file, nloop = 20)
    @test packed_sys isa Packmol.PackmolSystem
    @test packed_sys.status == :packing_ready
    @test isfile(tmp_output_file)
    rm(tmp_output_file, force=true)

end

@testitem "SolutionBoxUS tilde paths" begin
    using Packmol
    test_dir = Packmol.RecipesDirectory*"/test"
    if startswith(test_dir, homedir())
        tilde_dir = "~" * test_dir[length(homedir())+1:end]
        system = SolutionBoxUS(
            solute_pdbfile = "$tilde_dir/data/poly_h.pdb",
            solvent_pdbfile = "$tilde_dir/data/water.pdb",
            density=1.0,
        )
        @test system.solute_pdbfile == "$test_dir/data/poly_h.pdb"
        @test system.solvent_pdbfile == "$test_dir/data/water.pdb"
        input = tempname() * ".inp"
        redirect_stdout(devnull) do
            write_packmol_input(system; input, margin=5.0)
        end
        @test !any(contains("~"), readlines(input))
        @test Packmol.read_packmol_input(input).structure_types[1].filename == "$test_dir/data/poly_h.pdb"
    end
end

@testitem "Recipe keyword-only constructors" begin
    using Packmol
    # positional calls and missing mandatory keywords give an ArgumentError
    # listing the recipe's keywords, not a bare MethodError/UndefKeywordError
    for R in (SolutionBoxUS, SolutionBoxUSC, SolutionBoxUWI, Membrane, MembraneProtein)
        err = try R("protein.pdb") catch e e end
        @test err isa ArgumentError
        @test contains(err.msg, "keyword arguments only")
        @test contains(err.msg, "(; $(first(Packmol._mandatory_keywords(R)))=\"protein.pdb\", ...)")
        err = try R() catch e e end
        @test err isa ArgumentError
        @test contains(err.msg, "`$(first(Packmol._mandatory_keywords(R)))` is mandatory")
        @test contains(err.msg, "Optional keyword arguments:")
    end
end

@testitem "dodecahedral and octahedral cells: the solute within the margin from the faces" begin
    using PDBToolsBase: read_pdb
    using StaticArrays, LinearAlgebra
    atoms = read_pdb(joinpath(Packmol.RecipesDirectory, "test", "data", "poly_h.pdb"))
    x = [SVector(a.x, a.y, a.z) for a in atoms]
    x .-= Ref(sum(x) / length(x))
    for pbc in (:dodecahedral, :octahedral), margin in (5.0, 20.0)
        d = Packmol._wigner_seitz_size(atoms, margin, pbc)
        cell = pbc == :dodecahedral ? dodecahedral_unitcell(Float64, d) : octahedral_unitcell(Float64, d)
        # every atom at least the margin from every face: x⋅v̂ + margin ≤ |v|/2
        worst = maximum(dot(xi, v / norm(v)) + margin - norm(v) / 2
            for xi in x, v in (cell * SVector(i, j, k) for i in -2:2, j in -2:2, k in -2:2 if (i, j, k) != (0, 0, 0)))
        @test worst <= 1e-8
        @test worst >= -1e-6 # the smallest such cell: some atom touches the margin
    end
end
