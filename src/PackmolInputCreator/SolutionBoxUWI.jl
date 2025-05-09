mutable struct SolutionBoxUWI <: SolutionBox
    solute_pdbfile::String
    solute_charge::Int
    water_pdbfile::String
    cation_pdbfile::String
    anion_pdbfile::String
    cation_charge::Int
    anion_charge::Int
    solute_molar_mass::typeof(1.0u"g/mol")
    cation_molar_mass::typeof(1.0u"g/mol")
    anion_molar_mass::typeof(1.0u"g/mol")
    density_table::DensityTable
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
    solute_molar_mass::Union{Nothing,Number}=nothing,
    cation_molar_mass::Union{Nothing,Number}=nothing,
    anion_molar_mass::Union{Nothing,Number}=nothing,
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
        isnothing(cation_charge) && (cation_charge = 1)
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
    cation_molar_mass = cation_molar_mass * 1u"g/mol"

    if isnothing(anion_pdbfile)
        anion = [ Atom(index=1, name="CLA", resname="CLA", resnum=1, x=0.0f0, y=0.0f0, z=0.0f0) ]
        anion_pdbfile = scratch_dir * "/CLA.pdb"
        isnothing(anion_charge) && (anion_charge = -1)
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
    anion_molar_mass = anion_molar_mass * 1u"g/mol"

    solute = read_pdb(solute_pdbfile)
    if isnothing(solute_molar_mass)
        solute_molar_mass = mass(solute) * 1u"g/mol"
    end
    if isnothing(solute_charge)
        solute_charge = sum(charge, eachresidue(solute)) 
    end

    return SolutionBoxUWI(
        solute_pdbfile,
        solute_charge,
        water_pdbfile,
        cation_pdbfile,
        anion_pdbfile,
        cation_charge,
        anion_charge,
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
        Solute molar mass:  $(system.solute_molar_mass)
        Cation molar mass: $(system.cation_molar_mass)
        Anion molar mass: $(system.anion_molar_mass)
    ==================================================================
    """))
end

# Density of NaCl(aq) (concentrations in mol/kg, densities in g/mL), at 25°C
# From: https://advancedthermo.com/electrolytes/density_NaCl_Jun2021.html
const density_NaCl_aq = [ 
   0.1u"mol/kg" 	   1.00116u"g/mL" 
   0.2u"mol/kg" 	   1.00520u"g/mL" 
   0.3u"mol/kg" 	   1.00921u"g/mL" 
   0.4u"mol/kg" 	   1.01317u"g/mL" 
   0.5u"mol/kg" 	   1.01709u"g/mL" 
   0.6u"mol/kg" 	   1.02098u"g/mL" 
   0.7u"mol/kg" 	   1.02483u"g/mL" 
   0.8u"mol/kg" 	   1.02866u"g/mL" 
   0.9u"mol/kg" 	   1.03245u"g/mL" 
   1.0u"mol/kg" 	   1.03621u"g/mL" 
   1.2u"mol/kg" 	   1.04366u"g/mL" 
   1.4u"mol/kg" 	   1.05096u"g/mL" 
   1.6u"mol/kg" 	   1.05817u"g/mL" 
   1.8u"mol/kg" 	   1.06527u"g/mL" 
   2.0u"mol/kg" 	   1.07227u"g/mL" 
   2.5u"mol/kg" 	   1.08932u"g/mL" 
   3.0u"mol/kg" 	   1.10579u"g/mL" 
   3.5u"mol/kg" 	   1.12170u"g/mL" 
   4.0u"mol/kg" 	   1.13709u"g/mL" 
   4.5u"mol/kg" 	   1.15199u"g/mL" 
   5.0u"mol/kg" 	   1.16644u"g/mL" 
   5.5u"mol/kg" 	   1.18048u"g/mL" 
   6.0u"mol/kg" 	   1.19412u"g/mL" 
]

"""
    write_packmol_input(
        system::SolutionBoxUWI;
        ionic_concentration::Real=0.16u"mol/L",
        # Set density or provide a density table
        density::Union{Nothing,Number}=nothing, # g/mL
        density_table::Union{Nothing,AbstractMatrix{<:Number}}=density_NaCl_aq,
        # Optional set names
        input="box.inp",
        output="system.pdb",
        # box size or margin: default is a 20u"Å" margin
        box_sides::AbstractVector{<:Number}, # or
        margin::Number=20u"Å",
        # Optional: defaults to cubic box
        cubic::Bool=true,
    )

"""
function write_packmol_input(
    system::SolutionBoxUWI;
    ionic_concentration::Real=0.16, # mol/L
    input="box.inp",
    output="system.pdb",
    # box size
    box_sides::Union{Nothing,AbstractVector{<:Real}}=nothing, # or
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
    nions = round(Int, CMC * ionic_concentration * vs)

    # The actual number of cations or anions depends on their charges
    ncation = nions ÷ system.cation_charge - round(Int, system.solute_charge/2, RoundUp)
    nanion = nions ÷ (-1*system.anion_charge) + round(Int, system.solute_charge/2, RoundDown)

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
   
       Target ionic concentration = $ionic_concentration mol/L
   
       Box volume = $vbox Å³
       Solution volume (excluding solute) = $vs Å³   
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

    return
   
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


