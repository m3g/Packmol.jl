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
        Density units: $(system.density_units) ($(unit_name(system.density_units)))
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
    Mp = solute_molar_mass
    Mw = solvent_molar_mass

    # Density of pure solvent (g/mL)
    ρs = system.density
    
    # Molarity of pure solvent (mol/L)
    ms = 1000 * ρs / Mw

    # Convert solvent concentration in molecules/Å³
    cs = CMC * ms

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
    # as if it had the same mass density of the pure solvent
    vs = vbox - CMV * Mp / ρs

    # number of solvent molecules (molecules/Å³ * Å³)
    ns = round(Int, cs * vs)

    # Number of solvent molecules
    if ns == 0
        throw(ArgumentError("Number of solvent molecules is zero."))
    end

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $ms mol/L
        (of solvent   )      = $cs molecules/Å³
                             = $(system.density) g/mL

        Box volume = $vbox Å³
        Solution volume = $vs Å³   
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 

        Solute molar mass = $Mp g/mol
        Solvent molar mass = $Mw g/mol

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
            pbc $(join( -1.0*l, " ")), $(join( l, " "))

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
