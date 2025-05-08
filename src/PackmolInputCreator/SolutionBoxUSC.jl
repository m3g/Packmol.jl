struct DensityTable
    concentration_units::String
    concentration::Vector{<:Number}
    density::Vector{<:Number}
end

mutable struct SolutionBoxUSC <: SolutionBox
    solute_pdbfile::String
    solvent_pdbfile::String
    cossolvent_pdbfile::String
    density_table::DensityTable
    solute_molar_mass::Number
    solvent_molar_mass::Number
    cossolvent_molar_mass::Number
end

function DensityTable(
    concentration_units::String,
    concentration::AbstractVector{<:Number},
    density::AbstractVector{<:Number},
    M_solvent::Quantity,
    M_cossolvent::Quantity,
    rho_pure_solvent,
)
    cvec = typeof(1.0u"mol/L")[
        cconvert(c, concentration_units => "mol/L";
            M_solvent, M_solute=M_cossolvent, rho_solution=rho_pure_solvent * 1u"g/ml",
        ) for c in concentration 
    ]
    dvec = typeof(1.0u"g/mL")[uconvert(u"g/mL", 1u"g/mL" * d) for d in density]
    return DensityTable(concentration_units, cvec, dvec)
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
    if !(ustrip(density_table[begin, 1]) == 0.0)
        throw(ArgumentError("First line of density table must be the density of pure solvent, with cossolvent concentration equal to 0.0"))
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
        density_table[begin, 2],
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
        Density of pure solvent: $(first(system.density_table.density))
        Cosolvent concentration at highest concentration: $(last(system.density_table.concentration[end]))
        Density of cosolvent solution at highest concentration: $(last(system.density_table.density[end]))
        Molar masses: 
            solute: $(system.solute_molar_mass)
            solvent: $(system.solvent_molar_mass)
            cossolvent: $(system.cossolvent_molar_mass)
        Cocentration range: $(first(system.density_table.concentration)) - $(last(system.density_table.concentration))
    ==================================================================
    """))
end

"""
    write_packmol_input(
        system::SolutionBoxUSC;
        concentration::Real, 
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real,
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
    concentration::Real, 
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Real},Nothing} = nothing,
    margin::Union{Real,Nothing} = nothing, 
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

    # molar masses (g/mol)
    Mp = solute_molar_mass
    Mc = cossolvent_molar_mass
    Mw = solvent_molar_mass

    # Check consistency of the concentrations given
    cunit = system.concentration_units
    ρs = @view(system.density_table[:, 1])
    if cunit == "mol/L" && (first(ρs) < 1e-3 && last(ρs) ≈ 1)
        throw(ArgumentError("Concentrations in density table do not appear to be in mol/L."))
    end
    if (cunit in ("x", "mm", "vv")) && (any(x -> !(0 <= x <= 1), ρs)) 
        throw(ArgumentError("Concentrations in density table outside [0,1] range, and units are: $cunit"))
    end

    # Find the density corresponding to the target concentration
    ρ = interpolate_concentration(system, concentration)

    # Obtain the concentration in all units, for testing
    c_x = convert_concentration(system, concentration, cunit => "x")
    c_vv = convert_concentration(system, concentration, cunit => "vv")
    cc_mol = convert_concentration(system, concentration, cunit => "mol/L")
    cc_mm = convert_concentration(system, concentration, cunit => "mm")

    # aliases for clearer formulas
    ρc = density_pure_cossolvent(system) 
    ρw = density_pure_solvent(system) 

    # Convert cossolvent concentration in molecules/Å³
    cc = CMC * cc_mol

    # Set box side
    if isnothing(box_sides) && isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided."))
    elseif !isnothing(box_sides) && !isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided, but not both."))
    end
    solute_atoms = read_pdb(system.solute_pdbfile)
    solute_extrema = round.(maxmin(solute_atoms).xlength; digits=3)
    if !isnothing(margin)
        box_sides = solute_extrema .+ 2 .* margin
    end

    # Box volume (Å³)
    if !cubic
        vbox = box_sides[1] * box_sides[2] * box_sides[3]
    else
        max_side = maximum(box_sides)
        vbox = max_side^3
        box_sides .= max_side
    end

    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the solution
    vs = vbox - CMV * Mp / ρ

    # number of cossolvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, cc * vs)

    # Number of solvent molecules
    if nc != 0
        nw = round(Int, nc * (1 - c_x) / c_x)
    else
        nw = round(Int, vs * (ρw/CMV) / Mw)
    end

    # Final density of the solution (not inclusing solute volume)
    ρ = CMV * (Mc * nc + Mw * nw) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = 1000 * (nc / vs) * CMV

    # Final solvent concentration (mol/L)
    cw_f = 1000 * (nw / vs) * CMV

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $cc_mol mol/L
        (of cossolvent)      = $c_x molar fraction
                             = $c_vv volume fraction
                             = $cc_mm mass fraction 
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 

        Density of solution = $ρ g/mL
        Solute molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol
        Solvent molar mass = $Mw g/mol

        Number of cossolvent ($(basename(cossolvent_pdbfile))) molecules = $nc 
        Number of solvent ($(basename(solvent_pdbfile))) molecules = $nw 

        Final cossolvent concentration = $cc_f mol/L
                                    = $(CMC*cc_f) molecules/Å³
        Final solvent concentration = $cw_f mol/L
                                    = $(CMC*cw_f) molecules/Å³
                                    
        Final volume fraction = $((nc * Mc / ρc)/((nc * Mc / ρc) + (nw * Mw / ρw)))
        Final molar fraction = $(nc/(nc+nw))

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
            pbc $(join( -1.0*l, " ")), $(join( l, " "))

            structure $solute_pdbfile
                number 1
                center
                fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $solvent_pdbfile
                number $nw
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
        return nw, nc, 2*l
    else
        return nothing
    end
end # function write_packmol_input

@testitem "SolutionBoxUSC" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

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
