mutable struct SolutionBoxUSC <: SolutionBox
    solute_pdbfile::String
    solvent_pdbfile::String
    cossolvent_pdbfile::String
    density_table::DensityTable
    solute_molar_mass::Number
    solvent_molar_mass::Number
    cossolvent_molar_mass::Number
end

"""
    SolutionBoxUSC(; 
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        density_table::Matrix{<:Number},
        concentration_unit = "x",
        solute_molar_mass = nothing, # optional
        solvent_molar_mass = nothing, # optional
        cossolvent_molar_mass = nothing, # optional
    )

Setup a system composed of a solute (U) a solvent (S) and a cossolvent (C). 

The concentration units of the density table can be provided explicitly and
are assumed by default to be the molar fraction, `x`, of the cossolvent. The
density units are assumed to be `g/mL` unless the values are explicitly provided
with Unitful units. 

The molar masses of the solute, solvent, and cossolvent can be provided manually. If
not, they will be computed from the atom types in the PDB file.

"""
function SolutionBoxUSC(;
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        density_table::Matrix{<:Number},
        concentration_units::Union{Nothing,String} = nothing,
        solute_molar_mass::Union{Nothing,<:Number} = nothing,
        solvent_molar_mass::Union{Nothing,<:Number} = nothing,
        cossolvent_molar_mass::Union{Nothing,<:Number} = nothing,
    )
    if isnothing(concentration_units) && unit(density_table[begin, 1]) == NoUnits
        concentration_units = "x"
        @warn "Concentration units not provided, assuming molar fraction (x)." _file=nothing _line=nothing
    end
    if unit(density_table[begin, 2]) == NoUnits
        @warn "Density units not provided, assuming g/mL." _file=nothing _line=nothing
    end
    isnothing(solute_molar_mass) && (solute_molar_mass = mass(read_pdb(solute_pdbfile)))
    solute_molar_mass = _ensure_unit(solute_molar_mass, u"g/mol")
    isnothing(solvent_molar_mass) && (solvent_molar_mass = mass(read_pdb(solvent_pdbfile)))
    solvent_molar_mass = _ensure_unit(solvent_molar_mass, u"g/mol")
    isnothing(cossolvent_molar_mass) && (cossolvent_molar_mass = mass(read_pdb(cossolvent_pdbfile)))
    cossolvent_molar_mass = _ensure_unit(cossolvent_molar_mass, u"g/mol")
    # Convert concenrations in density table to mol/L
    dtable = DensityTable(
        concentration_units, 
        density_table[:, 1], 
        density_table[:, 2],
        solvent_molar_mass,
        cossolvent_molar_mass,
    )
    # Convert density table to g/mL
    system = SolutionBoxUSC(
        solute_pdbfile,
        solvent_pdbfile,
        cossolvent_pdbfile,
        dtable,
        solute_molar_mass,
        solvent_molar_mass,
        cossolvent_molar_mass,
    )
    return system
end

function Base.show(io::IO, ::MIME"text/plain", system::SolutionBoxUSC)
    print(io, chomp("""
    ==================================================================
    SolutionBoxUSC properties (Solute + Solvent + Cossolvent):
    ==================================================================
        Solute pdb file: $(basename(system.solute_pdbfile))
        Solvent pdb file: $(basename(system.solvent_pdbfile))
        Cossolvent pdb file: $(basename(system.cossolvent_pdbfile))
        Molar masses: 
            solute: $(system.solute_molar_mass)
            solvent: $(system.solvent_molar_mass)
            cossolvent: $(system.cossolvent_molar_mass)
        Cocentration range: $(first(system.density_table.concentration)) - $(last(system.density_table.concentration))
        Density range: $(first(system.density_table.density)) - $(last(system.density_table.density))
    ==================================================================
    """))
end

"""
    write_packmol_input(
        system::SolutionBoxUSC;
        concentration::Number, 
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        cubic::Bool = false,
    )

Function that generates an input file for Packmol. 

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternativelly, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

If `cubic` is set to true, the box will be cubic, and the box sides will be
equal in all 3 dimensions, respecting the minimum margin provided.

"""
function write_packmol_input(
    system::SolutionBoxUSC;
    concentration::Number, 
    concentration_units::Union{Nothing,String} = nothing,
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing} = nothing,
    margin::Union{Number,Nothing} = nothing, 
    cubic::Bool = false,
    # testing option
    debug = false,
)

    (; solute_pdbfile,
       solvent_pdbfile,
       cossolvent_pdbfile,
       solute_molar_mass,
       solvent_molar_mass,
       cossolvent_molar_mass,
    ) = system

    # molar masses
    Mu = solute_molar_mass
    Mc = cossolvent_molar_mass
    Ms = solvent_molar_mass

    # Check the input units of the concentration
    if unit(concentration) == NoUnits
        isnothing(concentration_units) && @warn "Concentration units not provided, assuming molar fraction (x)." _file=nothing _line=nothing
        cunit = "x"
    else
        cunit = UNIT_TYPE_MAP[unit(concentration)]
    end

    # Find the density corresponding to the target concentration
    ρ = interpolate_density(system, concentration, cunit)

    # Obtain the concentration in all units, for testing
    c_x = cconvert(concentration, cunit => "x"; M_solvent = Ms, M_solute = Mc, rho_solution = ρ)
    cc_mol = cconvert(concentration, cunit => "mol/L";  M_solvent = Ms, M_solute = Mc, rho_solution = ρ)
    cc_mm = cconvert(concentration, cunit => "w/w";  M_solvent = Ms, M_solute = Mc, rho_solution = ρ)

    # aliases for clearer formulas
    ρs = density_pure_solvent(system) 

    box_sides, solute_extrema = set_box_sides(system, box_sides, margin, cubic)
    vbox = prod(box_sides) 

    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the solution
    vs = vbox - uconvert(u"Å^3", Mu / ρs / Unitful.Na)

    # number of cossoUnitful.Na * concentration * vslvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, Unitful.Na * cc_mol * vs)

    # Number of solvent molecules
    if nc != 0
        ns = round(Int, nc * (1 - c_x) / c_x)
    else
        ns = round(Int, Unitful.Na * vs * ρs / Ms)
    end

    # Final density of the solution (not including solute volume)
    ρ = (Mc * nc + Ms * ns) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = uconvert(u"mol/L", (nc / Unitful.Na) / vs)

    # Final solvent concentration (mol/L)
    cs_f = uconvert(u"mol/L", (ns / Unitful.Na) / vs)

    # Half of box sides, to center the solute at the origin
    l = round.(typeof(1.0u"Å"), box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $cc_mol
        (of cossolvent)      = $c_x (molar fraction)
                             = $cc_mm (mass fraction)

        Box volume = $vbox
        Solution volume = $vs
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ]
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ]

        Density of solution = $ρ
        Solute molar mass = $Mu
        Cossolvent molar mass = $Mc
        Solvent molar mass = $Ms

        Number of cossolvent ($(basename(cossolvent_pdbfile))) molecules = $nc 
        Number of solvent ($(basename(solvent_pdbfile))) molecules = $ns 

        Final cossolvent concentration = $cc_f 
                                    = $(inv(cconvert(cc_f, "mol/L" => "NumberDensity"))) per molecule.
        Final solvent concentration = $cs_f
                                    = $(inv(cconvert(cs_f, "mol/L" => "NumberDensity"))) per molecule.
                                    
        Final molar fraction = $(nc/(nc+ns))

        Cubic box requested: $cubic

        ==================================================================
        """
    println(summary)

    open(input, "w") do io
        print(io,
            """
            # 
            # Packmol input file
            # 
            # Generated by MolSimToolkit.jl
            #
            """
        )
        for line in split(summary, "\n")
            println(io, "# $line")
        end
        println(io,
            """
            #
            tolerance 2.0
            output $output
            add_box_sides 1.0
            filetype pdb
            seed -1
            packall
            pbc $(join(-1.0*ustrip(l), " ")), $(join(ustrip(l), " "))

            structure $solute_pdbfile
                number 1
                center
                fixed 0.0 0.0 0.0 0.0 0.0 0.0
            end structure

            structure $solvent_pdbfile
                number $ns
            end structure
            """)
        if nc > 0
            println(io,
                """
                structure $cossolvent_pdbfile
                    number $nc
                end structure
                """)
        end
    end
    print(chomp(
        """
        Wrote file: $input

        ==================================================================
        """))
    
    if debug 
        return ns, nc, 2*l
    else
        return nothing
    end
end # function write_packmol_input

@testitem "SolutionBoxUSC" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.PackmolInputCreatorDirectory*"/test"

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
            solute: 5612.791939999981 g mol^-1
            solvent: 46.06922 g mol^-1
            cossolvent: 18.01534 g mol^-1
        Cocentration range: 0.0 mol L^-1 - 55.402784515862606 mol L^-1
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
