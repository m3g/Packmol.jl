mutable struct SolutionBoxUWI <: SolutionBox
    density_units::String
    solute_pdbfile::String
    solute_charge::Int
    water_pdbfile::String
    cation_pdbfile::String
    anion_pdbfile::String
    cation_charge::Int
    anion_charge::Int
    density::Float64
    solute_molar_mass::Float64
    cation_molar_mass::Float64
    anion_molar_mass::Float64
end

"""
    SolutionBoxUWI(; 
        solute_pdbfile::String, 
        solute_charge::Integer, # optional - computed from PDB file by default
        water_pdbfile::String, # optional
        cation_pdbfile::String, # optional - SOD by default
        anion_pdbfile::String, # optional - CLA by default
        cation_charge::Integer, # optional - 1 by default
        anion_charge::Integer, # optional - -1 by default
        density::Real=1.00, # optional - 1.00 g/mL by default
        density_units::String = "g/mL", # optional
        solute_molar_mass = nothing, # optional
        cation_molar_mass = nothing, # optional
        anion_molar_mass = nothing, # optional
    )
    
Setup a system composed of a solute (U), water (W), cation (I) and anion (I).

"""
function SolutionBoxUWI(;
    solute_pdbfile::String,
    solute_charge::Union{Nothing,Integer}=nothing,
    water_pdbfile::Union{Nothing,String}=nothing,
    cation_pdbfile::Union{Nothing,String}=nothing,
    anion_pdbfile::Union{Nothing,String}=nothing,
    cation_charge::Union{Nothing,Integer}=nothing,
    anion_charge::Union{Nothing,Integer}=nothing,
    density::Real=1.0,
    density_units::Union{Nothing,String}="g/mL",
    solute_molar_mass::Union{Nothing,Real}=nothing,
    cation_molar_mass::Union{Nothing,Real}=nothing,
    anion_molar_mass::Union{Nothing,Real}=nothing,
)
    scratch_dir = tempname()
    mkdir(scratch_dir)

    if isnothing(water_pdbfile)
        water = [
            Atom(index=1, name="OH2", resname="HOH", resnum=1, x=0.0f0, y=0.0f0, z=0.0f0),
            Atom(index=2, name="H1", resname="HOH", resnum=1, x=0.758602f0, y=0.0f0, z=0.504284f0),
            Atom(index=3, name="H2", resname="HOH", resnum=1, x=0.758602f0, y=0.0f0, z=-0.504284f0), 
        ]
        water_pdbfile = scratch_dir * "/HOH.pdb"
        write_pdb(water_pdbfile, water)
    end
    if isnothing(cation_pdbfile)
        cation = [ Atom(index=1, name="SOD", resname="SOD", resnum=1, x=0.0f0, y=0.0f0, z=0.0f0) ]
        cation_pdbfile = scratch_dir * "/SOD.pdb" 
        cation_charge = 1
        write_pdb(cation_pdbfile, cation)
    else
        cation = read_pdb(cation_pdbfile)
        if isnothing(cation_charge)
            throw(ArgumentError("""\n
                Custom cation PDB file provided, but no charge specified.
                Use `cation_charge` keyword argument to specify the charge of the cation.

            """))
        end
    end
    if isnothing(cation_molar_mass)
        cation_molar_mass = mass(cation)
    end

    if isnothing(anion_pdbfile)
        anion = [ Atom(index=1, name="CLA", resname="CLA", resnum=1, x=0.0f0, y=0.0f0, z=0.0f0) ]
        anion_pdbfile = scratch_dir * "/CLA.pdb"
        anion_charge = -1
        write_pdb(anion_pdbfile, anion)
    else
        anion = read_pdb(anion_pdbfile)
        if isnothing(anion_charge)
            throw(ArgumentError("""\n
                Custom anion PDB file provided, but no charge specified.
                Use `anion_charge` keyword argument to specify the charge of the anion.

            """))
        end
    end
    if isnothing(anion_molar_mass)
        anion_molar_mass = mass(anion)
    end

    solute = read_pdb(solute_pdbfile)
    if isnothing(solute_molar_mass)
        solute_molar_mass = mass(solute)
    end
    if isnothing(solute_charge)
        solute_charge = sum(charge, eachresidue(solute)) 
    end

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
    # Convert density in mol/L to g/mL
    if density_units == "mol/L"
        density = density * solvent_molar_mass / 1000
    end

    return SolutionBoxUWI(
        density_units,
        solute_pdbfile,
        solute_charge,
        water_pdbfile,
        cation_pdbfile,
        anion_pdbfile,
        cation_charge,
        anion_charge,
        density,
        solute_molar_mass,
        cation_molar_mass,
        anion_molar_mass,
    )
end

function Base.show(io::IO, ::MIME"text/plain", system::SolutionBoxUWI)
    print(io, chomp("""
    ==================================================================
    SolutionBoxUWI properties (Solute + Water + Ions):
    ==================================================================
        Solute pdb file: $(basename(system.solute_pdbfile))
        Solute charge: $(system.solute_charge)
        Water pdb file: $(basename(system.water_pdbfile))
        Cation pdb file and charge: $(basename(system.cation_pdbfile)) +$(system.cation_charge)
        Anion pdb file: $(basename(system.anion_pdbfile)) $(system.anion_charge)
        Density: $(system.density) g/mL
        Ion concentration: $(system.ionic_concentration) mol/L
        Solute molar mass:  $(system.solute_molar_mass) g/mol
        Cation molar mass: $(system.cation_molar_mass) g/mol
        Anion molar mass: $(system.anion_molar_mass) g/mol
    ==================================================================
    """))
end



"""
    write_packmol_input(
        system::SolutionBoxUSC;
        ionic_concentration::Real=0.16, # mol/L
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real,
        cubic::Bool = false,
    )

"""
function write_packmol_input(
    system::SolutionBoxUWI;
    ionic_concentration::Real=0.16, # mol/L
    input="box.inp",
    output="system.pdb",
    # box size
    box_sides::AbstractVector{<:Real}=nothing, # or
    margin::Real=0.0,
    cubic::Bool = false,
)
    if isnothing(box_sides) && isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided."))
    end
    if !isnothing(box_sides) && !isnothing(margin)
        throw(ArgumentError("Only one of box_sides or margin can be provided."))
    end
    if !isnothing(box_sides) && length(box_sides) != 3
        throw(ArgumentError("box_sides must be a vector of length 3."))
    end
    if !isnothing(margin) && margin <= 0.0
        throw(ArgumentError("Margin must be positive."))
    end

    # Set box side
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

    # Estimated solute volume
    vsolute = CMV * system.solute_molar_mass / system.density

    # Solution volume
    vs = vbox - vsolute

    # Number of ionic "molecules" in solution
    nions = round(Int, ionic_concentration * vs / CMC)

    # The actual number of cations or anions depends on their charges
    ncation = nions ÷ system.cation_charge - round(Int, system.solute_charge/2, RoundUp)
    nanion = nions ÷ system.anion_charge + round(Int, system.solute_charge/2, RoundDown)

    # Total final charge (for checking)
    total_charge = ncation * system.cation_charge + nanion * system.anion_charge + system.solute_charge

    # Total ion mass
    mions = ncation * system.cation_molar_mass + nanion * system.anion_molar_mass

    # Water mass
    mwater = system.density * vs - mions

    # Number of water molecules
    nwater = round(Int, mwater / system.cation_molar_mass)

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    summary = """
       ==================================================================
       Summary:
       ==================================================================
   
       Target ionic concentration = $(system.ionic_concentration) mol/L
   
       Box volume = $vbox Å³
       Solution volume = $vs Å³   
       Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
       Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 
   
       Number of cations = $ncation
       Number of anions = $nanion
       Solute charge = $(system.solute_charge)
       Total system charge = $total_charge
       Number water molecules = $nwater
   
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

end


