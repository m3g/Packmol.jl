mutable struct SolutionBoxUS <: SolutionBox
    solute_pdbfile::String
    solvent_pdbfile::String
    density::Quantity
    solute_molar_mass::Quantity
    solvent_molar_mass::Quantity
end

"""
    SolutionBoxUS(; 
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        density::Union{Quantity,Real}, # density, using Unitful units, or assumed to be g/mL
        solute_molar_mass = nothing, # optional
        solvent_molar_mass = nothing, # optional
    )

Setup a system composed of a solute (U) and a solvent (S). 

The mass or molar density can be provided. If unitless specified, the density will be assumed to be in g/mL, and the molar masses in g/mol.

If the molar masses are not provided, they will be computed from the atom types in the PDB file.

"""
function SolutionBoxUS(;
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        density::Number, # density, using Unitful units, or assumed to be g/mL
        solute_molar_mass::Union{Nothing,Number} = nothing,
        solvent_molar_mass::Union{Nothing,Number} = nothing,
    )
    if unit(density) == NoUnits
        @warn "Density units not provided, assuming g/mL." _file=nothing _line=nothing
        density = density * 1.0u"g/mL"
    end
    ustrip(density) <= 0.0 && throw(ArgumentError("Density must be positive."))
    isnothing(solute_molar_mass) && (solute_molar_mass = mass(read_pdb(solute_pdbfile)) * 1.0u"g/mol")
    solute_molar_mass = _ensure_unit(solute_molar_mass, u"g/mol")
    isnothing(solvent_molar_mass) && (solvent_molar_mass = mass(read_pdb(solvent_pdbfile)) * 1.0u"g/mol")
    solvent_molar_mass = _ensure_unit(solvent_molar_mass, u"g/mol")
    # Convert density in mol/L to g/mL
    unit(density) == u"mol/L" && (density = uconvert(u"g/mL", density * solvent_molar_mass))
    # Construct system
    system = SolutionBoxUS(
        solute_pdbfile,
        solvent_pdbfile,
        density,
        solute_molar_mass,
        solvent_molar_mass,
    )
    return system
end

function Base.show(io::IO, ::MIME"text/plain", system::SolutionBoxUS)
    print(io, chomp("""
    ==================================================================
    SolutionBoxUS properties (Solute + Solvent):
    ==================================================================
        Solute pdb file: $(basename(system.solute_pdbfile))
        Solvent pdb file: $(basename(system.solvent_pdbfile))
        Density of pure solvent: $(system.density)
        Molarity of pure solvent: $(uconvert(u"mol/L", system.density / system.solvent_molar_mass))
        Molar masses: 
            solute: $(system.solute_molar_mass)
            solvent: $(system.solvent_molar_mass)
    ==================================================================
    """))
end

"""
    write_packmol_input(
        system::SolutionBoxUS;
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        cubic::Bool = false,
    )

Function that generates an input file for Packmol for a Solute + Solvent system.

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternativelly, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

If `cubic` is set to true, the box will be cubic, and the box sides will be
equal in all 3 dimensions, respecting the minimum margin provided.

"""
function write_packmol_input(
    system::SolutionBoxUS;
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing} = nothing,
    margin::Union{<:Number,Nothing} = nothing, 
    cubic::Bool = false,
    # testing option
    debug = false,
)

    (; solute_pdbfile,
       solvent_pdbfile,
       solute_molar_mass,
       solvent_molar_mass,
    ) = system

    # molar masses (g/mol)
    Mp = solute_molar_mass
    Mw = solvent_molar_mass

    # Density of pure solvent (g/mL)
    ρs = system.density
    
    # Molarity of pure solvent (mol/L)
    ms = uconvert(u"mol/L", ρs / Mw)

    # Convert solvent concentration in molecules/Å³
    cs = cconvert(ms, "mol/L" => "molecules/Å^-3")

    # Set box side
    if isnothing(box_sides) && isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided."))
    elseif !isnothing(box_sides) && !isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided, but not both."))
    end
    solute_atoms = read_pdb(system.solute_pdbfile)
    solute_extrema = 1.0u"Å" * round.(maxmin(solute_atoms).xlength; digits=3)
    if !isnothing(margin)
        margin = _ensure_unit(margin, u"Å")
        box_sides = (solute_extrema .+ 2 .* margin)
    end
    box_sides = _ensure_unit.(box_sides, u"Å")

    # Box volume (Å³)
    if !cubic
        vbox = box_sides[1] * box_sides[2] * box_sides[3]
    else
        max_side = maximum(box_sides)
        vbox = max_side^3
        box_sides .= max_side
    end
        
    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the pure solvent
    vs = vbox - uconvert(u"Å^3", Mp / ρs / Unitful.Na)

    # number of solvent molecules (molecules/Å³ * Å³)
    ns = round(Int, cs * vs)

    # Number of solvent molecules
    ns == 0 && throw(ArgumentError("Number of solvent molecules is zero."))

    # Half of box sides, to center the solute at the origin
    l = round.(typeof(1.0u"Å"), box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $ms
        (of solvent   )      = $cs
                             = $ρs

        Box volume = $vbox
        Solution volume = $vs
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ]
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ]

        Solute molar mass = $Mp
        Solvent molar mass = $Mw 

        Number of solvent ($(basename(solvent_pdbfile))) molecules = $ns 

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
            pbc $(join( -1.0*ustrip(l), " ")), $(join(ustrip(l), " "))

            structure $solute_pdbfile
                number 1
                center
                fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $solvent_pdbfile
                number $ns
            end structure
            """)
    end
    print(chomp(
        """
        Wrote file: $input

        ==================================================================
        """))
    
    if debug 
        return ns, 2*l
    else
        return nothing
    end
end # function write_packmol_input

@testitem "SolutionBoxUS" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.PackmolInputCreatorDirectory*"/test"

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

end
