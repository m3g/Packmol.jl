mutable struct SolutionBoxUS <: SolutionBox
    density_units::String
    solute_pdbfile::String
    solvent_pdbfile::String
    density::Float64
    solute_molar_mass::Float64
    solvent_molar_mass::Float64
end

"""
    SolutionBoxUS(; 
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        density::Float64,
        density_units = "g/mL",
        solute_molar_mass = nothing, # optional
        solvent_molar_mass = nothing, # optional
    )

Setup a system composed of a solute (U) and a solvent (S). 

The mass or molar density can be provided, and the units are either
`"g/mL"` or `"mol/L"`. The default is `"g/mL"`. 

The molar masses of the solute and solvent can be provided manually. If
not, they will be computed from the atom types in the PDB file.

"""
function SolutionBoxUS(;
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        density::Float64,
        density_units::Union{Nothing,String} = nothing,
        solute_molar_mass::Union{Nothing,Real} = nothing,
        solvent_molar_mass::Union{Nothing,Real} = nothing,
    )
    if isnothing(density_units)
        density_units = "g/mL"
        @warn "Density units not provided, assuming g/mL." _file=nothing _line=nothing
    end
    if density <= 0.0
        throw(ArgumentError("Density must be positive."))
    end
    if !(density_units in ("g/mL", "mol/L"))
        throw(ArgumentError("Density units must be g/mL or mol/L."))
    end
    if isnothing(solute_molar_mass)
        solute_molar_mass = mass(read_pdb(solute_pdbfile))
    end
    if isnothing(solvent_molar_mass)
        solvent_molar_mass = mass(read_pdb(solvent_pdbfile))
    end
    # Convert density in mol/L to g/mL
    if density_units == "mol/L"
        density = density * solvent_molar_mass / 1000
    end
    system = SolutionBoxUS(
        density_units,
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
        Density of pure solvent: $(system.density) g/mL
        Molarity of pure solvent: $(1000 * system.density / system.solvent_molar_mass) mol/L
        Density units: $(system.density_units)
        Molar masses: 
            solute: $(system.solute_molar_mass) g/mol
            solvent: $(system.solvent_molar_mass) g/mol
    ==================================================================
    """))
end

"""
    write_packmol_input(
        system::SolutionBoxUS;
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real,
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
    box_sides::Union{AbstractVector{<:Real},Nothing} = nothing,
    margin::Union{Real,Nothing} = nothing, 
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
    Mp = solute_molar_mass * 1u"g/mol"
    Mw = solvent_molar_mass * 1u"g/mol"

    # Density of pure solvent (g/mL)
    ρs = system.density * 1u"g/mL"
    
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
        margin = margin * 1.0u"Å"
        box_sides = (solute_extrema .+ 2 .* margin)
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
    using ShowMethodTesting

    test_dir = Packmol.PackmolInputCreatorDirectory*"/test"

    # system with water only
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=55.5,
        density_units="mol/L",
    ) 
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.density ≈ 1.0 atol = 0.01
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
        density_units="g/mL",
    ) 
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.density ≈ 1.0 atol = 0.01
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
    ) 
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.density ≈ 1.0 atol = 0.01
    @test_throws ArgumentError SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=0.0,
    ) 
    @test_throws ArgumentError SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1000.0,
        density_units="g/L"
    ) 
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=55.5,
        density_units="mol/L",
        solvent_molar_mass=18.01534,
        solute_molar_mass=5612.79194,
    ) 
    @test system.density ≈ 1.0 atol = 0.01

    # Test show methods
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
        density_units="g/mL",
    ) 
    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUS properties (Solute + Solvent):
    ==================================================================
        Solute pdb file: poly_h.pdb
        Solvent pdb file: water.pdb
        Density of pure solvent: 1.0 g/mL
        Molarity of pure solvent: 55.508250191225926 mol/L
        Density units: g/mL
        Molar masses: 
            solute: 5612.791939999981 g/mol
            solvent: 18.01534 g/mol
    ==================================================================
    """

    # System with water only
    system = SolutionBoxUS(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        density=1.0,
        density_units="g/mL",
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

end
