# Reference density of pure water at 25°C, used only as a coarse volume
# estimate for the solute (mirroring how SolutionBoxUS/SolutionBoxUSC use the
# pure-solvent density for the same purpose) and as the density assumed when
# the target ionic concentration falls below the density table's range (i.e.
# indistinguishable from pure water for that purpose).
const PURE_WATER_DENSITY = 0.9970749u"g/mL"

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

mutable struct SolutionBoxUWI <: Recipe
    solute_pdbfile::String
    solute_charge::Int
    water_pdbfile::String
    cation_pdbfile::String
    anion_pdbfile::String
    cation_charge::Int
    anion_charge::Int
    solute_molar_mass::typeof(1.0u"g/mol")
    water_molar_mass::typeof(1.0u"g/mol")
    cation_molar_mass::typeof(1.0u"g/mol")
    anion_molar_mass::typeof(1.0u"g/mol")
    density_table::DensityTable
end

"""
    SolutionBoxUWI(;
        solute_pdbfile::String,
        solute_charge::Integer, # optional - computed from PDB file by default
        water_pdbfile::String, # optional
        cation_pdbfile::String, # optional - SOD (Na+) by default
        anion_pdbfile::String, # optional - CLA (Cl-) by default
        cation_charge::Integer, # optional - 1 by default
        anion_charge::Integer, # optional - -1 by default
        solute_molar_mass = nothing, # optional
        water_molar_mass = nothing, # optional
        cation_molar_mass = nothing, # optional
        anion_molar_mass = nothing, # optional
        density_table::Matrix{<:Number} = density_NaCl_aq, # optional
        density_table_units = "mol/kg", # optional
    )

Setup a system composed of a solute (U), water (W), a cation and an anion (I, for ions).

`cation_charge` must be a positive integer and `anion_charge` a negative integer. Molar
masses are computed from the atom types in the PDB files if not provided.

`density_table` provides the density of the salt solution as a function of ionic
concentration, used to size the number of water molecules. It defaults to the density of
aqueous NaCl, appropriate for the default (sodium/chloride) ions; a custom table should be
provided (as a `concentration_units => density` matrix, with `density_table_units` set
accordingly) if custom ions are used and an accurate solution density is required.

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
    water_molar_mass::Union{Nothing,Number}=nothing,
    cation_molar_mass::Union{Nothing,Number}=nothing,
    anion_molar_mass::Union{Nothing,Number}=nothing,
    density_table::Union{Nothing,AbstractMatrix{<:Number}}=nothing,
    density_table_units::String="mol/kg",
)
    custom_ions = !isnothing(cation_pdbfile) || !isnothing(anion_pdbfile)
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
    isnothing(water_molar_mass) && (water_molar_mass = mass(read_pdb(water_pdbfile)))
    water_molar_mass = _ensure_unit(water_molar_mass, u"g/mol")

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
    cation_charge <= 0 && throw(ArgumentError("cation_charge must be a positive integer, got $cation_charge."))
    isnothing(cation_molar_mass) && (cation_molar_mass = mass(cation))
    cation_molar_mass = _ensure_unit(cation_molar_mass, u"g/mol")

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
    anion_charge >= 0 && throw(ArgumentError("anion_charge must be a negative integer, got $anion_charge."))
    isnothing(anion_molar_mass) && (anion_molar_mass = mass(anion))
    anion_molar_mass = _ensure_unit(anion_molar_mass, u"g/mol")

    solute = read_pdb(solute_pdbfile)
    isnothing(solute_molar_mass) && (solute_molar_mass = mass(solute))
    solute_molar_mass = _ensure_unit(solute_molar_mass, u"g/mol")
    isnothing(solute_charge) && (solute_charge = round(Int, sum(charge, eachresidue(solute))))

    if isnothing(density_table)
        if custom_ions
            @warn "Custom ions provided without a `density_table`; assuming the density of aqueous NaCl. " *
                  "Provide `density_table` (and `density_table_units`) explicitly for an accurate density " *
                  "estimate with these ions." _file=nothing _line=nothing
        end
        density_table = density_NaCl_aq
    end
    dtable = DensityTable(
        density_table_units,
        density_table[:, 1],
        density_table[:, 2],
        water_molar_mass,
        cation_molar_mass + anion_molar_mass,
    )

    return SolutionBoxUWI(
        solute_pdbfile,
        solute_charge,
        water_pdbfile,
        cation_pdbfile,
        anion_pdbfile,
        cation_charge,
        anion_charge,
        solute_molar_mass,
        water_molar_mass,
        cation_molar_mass,
        anion_molar_mass,
        dtable,
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
        Anion pdb file and charge: $(basename(system.anion_pdbfile)) $(system.anion_charge)
        Molar masses:
            solute: $(system.solute_molar_mass)
            water: $(system.water_molar_mass)
            cation: $(system.cation_molar_mass)
            anion: $(system.anion_molar_mass)
        Ionic concentration range: $(first(system.density_table.concentration)) - $(last(system.density_table.concentration))
        Density range: $(first(system.density_table.density)) - $(last(system.density_table.density))
    ==================================================================
    """))
end

"""
    write_packmol_input(
        system::SolutionBoxUWI;
        ionic_concentration::Number=0.15u"mol/L",
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        cubic::Bool = false,
    )

Function that generates an input file for Packmol for a Solute + Water + Ions system.

`ionic_concentration` is the target bulk salt concentration (e.g. NaCl), assumed to be in
`mol/L` if given unitless. Enough additional counter-ions (of a single sign) are added, on
top of the bulk salt, to exactly neutralize the solute's own charge (`system.solute_charge`);
an `ArgumentError` is thrown if that charge cannot be exactly neutralized with the ions'
charges (e.g. neutralizing an odd solute charge with only doubly-charged counter-ions).

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternatively, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

If `cubic` is set to true, the box will be cubic, and the box sides will be
equal in all 3 dimensions, respecting the minimum margin provided.

"""
#
# Shared box-sizing/molecule-count/charge-neutralization computation for
# SolutionBoxUWI, used by both `write_packmol_input` (which prints/writes
# `summary` to a `.inp` file) and `packmol` (which builds a `PackmolSystem`
# directly and never writes one).
#
function _setup(
    system::SolutionBoxUWI,
    ionic_concentration::Number,
    box_sides::Union{AbstractVector{<:Number},Nothing},
    margin::Union{<:Number,Nothing},
    cubic::Bool,
)
    (; solute_charge,
       water_pdbfile,
       cation_pdbfile,
       anion_pdbfile,
       cation_charge,
       anion_charge,
       solute_molar_mass,
       water_molar_mass,
       cation_molar_mass,
       anion_molar_mass,
       density_table,
    ) = system

    if unit(ionic_concentration) == NoUnits
        @warn "Ionic concentration units not provided, assuming mol/L." _file=nothing _line=nothing
        ionic_concentration = ionic_concentration * 1.0u"mol/L"
    end
    ionic_concentration = _ensure_unit(ionic_concentration, u"mol/L")
    ustrip(ionic_concentration) < 0 && throw(ArgumentError("Ionic concentration must be non-negative."))

    # Molar mass of the (neutral) salt formula unit, e.g. NaCl
    Msalt = cation_molar_mass + anion_molar_mass

    # Set box sides and volume (reuses the same logic as SolutionBoxUS/SolutionBoxUSC)
    box_sides, solute_extrema = set_box_sides(system, box_sides, margin, cubic)
    vbox = prod(box_sides)

    # Solution volume (vbox - vsolute) - vsolute is estimated as if the solute
    # had the density of pure water
    vsolute = uconvert(u"Å^3", solute_molar_mass / PURE_WATER_DENSITY / Unitful.Na)
    vs = vbox - vsolute
    vs <= zero(vs) && throw(ArgumentError("Box is too small to fit the solute (increase margin or box_sides)."))

    # Number of bulk salt formula units for the target ionic concentration
    nsalt = round(Int, Unitful.Na * ionic_concentration * vs)

    # Charge-neutral bulk salt stoichiometry (handles asymmetric valences, e.g. CaCl₂:
    # 1 cation of charge +2 needs 2 anions of charge -1 to remain neutral on its own)
    ncation = nsalt * abs(anion_charge)
    nanion = nsalt * abs(cation_charge)

    # Extra counter-ions (of a single sign) to exactly neutralize the solute's own charge
    if solute_charge > 0
        if mod(solute_charge, abs(anion_charge)) != 0
            throw(ArgumentError(
                "Solute charge ($solute_charge) is not a multiple of the anion charge magnitude " *
                "($(abs(anion_charge))); cannot exactly neutralize the system with these ions."
            ))
        end
        nanion += solute_charge ÷ abs(anion_charge)
    elseif solute_charge < 0
        if mod(-solute_charge, cation_charge) != 0
            throw(ArgumentError(
                "Solute charge ($solute_charge) is not a multiple of the cation charge " *
                "($cation_charge); cannot exactly neutralize the system with these ions."
            ))
        end
        ncation += (-solute_charge) ÷ cation_charge
    end

    total_charge = ncation * cation_charge + nanion * anion_charge + solute_charge
    @assert total_charge == 0 "internal error: system is not charge-neutral (please report this)"

    # Solution density at the target ionic concentration (extra neutralizing ions are
    # assumed too dilute to shift the density noticeably). Below the density table's
    # range, the solution is indistinguishable from pure water for this purpose.
    ρ = if ionic_concentration <= first(density_table.concentration)
        PURE_WATER_DENSITY
    else
        interpolate_density(density_table, water_molar_mass, Msalt, ionic_concentration, "mol/L")
    end

    # Ion mass, and water mass filling the remaining solution mass
    mions = (ncation / Unitful.Na) * cation_molar_mass + (nanion / Unitful.Na) * anion_molar_mass
    msolution = uconvert(u"g", ρ * vs)
    mwater = msolution - uconvert(u"g", mions)
    mwater <= zero(mwater) && throw(ArgumentError(
        "Ion mass exceeds the solution mass: reduce ionic_concentration or enlarge the box."
    ))
    nwater = round(Int, Unitful.Na * mwater / water_molar_mass)

    # Half of box sides, to center the solute at the origin
    l = round.(typeof(1.0u"Å"), box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target ionic concentration = $ionic_concentration
        Solution density (at that concentration) = $ρ

        Box volume = $vbox
        Solution volume = $vs
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ]
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ]

        Solute molar mass = $solute_molar_mass
        Solute charge = $solute_charge

        Number of water ($(basename(water_pdbfile))) molecules = $nwater
        Number of cations ($(basename(cation_pdbfile)), charge +$cation_charge) = $ncation
        Number of anions ($(basename(anion_pdbfile)), charge $anion_charge) = $nanion
        Total system charge (check) = $total_charge

        Cubic box requested: $cubic

        ==================================================================
        """
    return (; nwater, ncation, nanion, l, summary)
end

"""
    write_packmol_input(
        system::SolutionBoxUWI;
        ionic_concentration::Number=0.15u"mol/L",
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        cubic::Bool = false,
    )

Function that generates an input file for Packmol for a Solute + Water + Ions system.

`ionic_concentration` is the target bulk salt concentration (e.g. NaCl), assumed to be in
`mol/L` if given unitless. Enough additional counter-ions (of a single sign) are added, on
top of the bulk salt, to exactly neutralize the solute's own charge (`system.solute_charge`);
an `ArgumentError` is thrown if that charge cannot be exactly neutralized with the ions'
charges (e.g. neutralizing an odd solute charge with only doubly-charged counter-ions).

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternatively, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

If `cubic` is set to true, the box will be cubic, and the box sides will be
equal in all 3 dimensions, respecting the minimum margin provided.

"""
function write_packmol_input(
    system::SolutionBoxUWI;
    ionic_concentration::Number=0.15,
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing}=nothing,
    margin::Union{<:Number,Nothing}=nothing,
    cubic::Bool=false,
    # testing option
    debug=false,
)
    (; solute_pdbfile, water_pdbfile, cation_pdbfile, anion_pdbfile) = system
    (; nwater, ncation, nanion, l, summary) = _setup(system, ionic_concentration, box_sides, margin, cubic)
    println(summary)

    open(input, "w") do io
        print(io,
            """
            #
            # Packmol input file
            #
            # Generated by Packmol.jl
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
            pbc $(join( -1.0*ustrip(l), " ")) $(join(ustrip(l), " "))

            structure $solute_pdbfile
                number 1
                center
                fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $water_pdbfile
                number $nwater
            end structure
            """)
        if ncation > 0
            println(io,
                """
                structure $cation_pdbfile
                    number $ncation
                end structure
                """)
        end
        if nanion > 0
            println(io,
                """
                structure $anion_pdbfile
                    number $nanion
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
        return nwater, ncation, nanion, 2*l
    else
        return nothing
    end
end # function write_packmol_input

"""
    packmol(
        system::SolutionBoxUWI;
        ionic_concentration::Number=0.15u"mol/L",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        cubic::Bool = false,
        kwargs...,
    )

Builds and packs a Solute + Water + Ions system directly, entirely in memory: equivalent
to calling [`write_packmol_input`](@ref write_packmol_input(::SolutionBoxUWI)) followed by
`packmol` on the resulting file, except no `.inp` file is ever written.

`ionic_concentration`, `output`, `box_sides`, `margin`, and `cubic` behave as in
`write_packmol_input`. Any other keyword (`nloop`, `iprint`, `seed`, `optimizer`, ...) is
forwarded to the packing engine — see `packmol(::PackmolSystem)`.

"""
function packmol(
    system::SolutionBoxUWI;
    ionic_concentration::Number=0.15,
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing}=nothing,
    margin::Union{<:Number,Nothing}=nothing,
    cubic::Bool=false,
    kwargs...,
)
    (; solute_pdbfile, water_pdbfile, cation_pdbfile, anion_pdbfile) = system
    (; nwater, ncation, nanion, l) = _setup(system, ionic_concentration, box_sides, margin, cubic)
    structure_types = [
        _fixed_solute_structure_type(solute_pdbfile),
        structure_type(water_pdbfile; number=nwater),
    ]
    ncation > 0 && push!(structure_types, structure_type(cation_pdbfile; number=ncation))
    nanion > 0 && push!(structure_types, structure_type(anion_pdbfile; number=nanion))
    packmol_system = PackmolSystem(structure_types;
        output, tolerance=2.0, add_box_sides=true, seed=-1, _recipe_unitcell(l)...,
    )
    return packmol(packmol_system; kwargs...)
end

@testitem "SolutionBoxUWI" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    test_dir = Packmol.RecipesDirectory*"/test"

    # Default ions (SOD/CLA), neutral solute
    system = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 0)
    @test system.cation_charge == 1
    @test system.anion_charge == -1
    @test system.cation_molar_mass ≈ 22.99u"g/mol" rtol = 0.01
    @test system.anion_molar_mass ≈ 35.45u"g/mol" rtol = 0.01
    @test system.water_molar_mass ≈ 18.02u"g/mol" rtol = 0.01

    @test parse_show(system) ≈ """
    ==================================================================
    SolutionBoxUWI properties (Solute + Water + Ions):
    ==================================================================
        Solute pdb file: poly_h.pdb
        Solute charge: 0
        Water pdb file: HOH.pdb
        Cation pdb file and charge: SOD.pdb +1
        Anion pdb file and charge: CLA.pdb -1
        Molar masses:
            solute: 5612.80078125 g mol^-1
            water: 18.01534080505371 g mol^-1
            cation: 22.989770889282227 g mol^-1
            anion: 35.452999114990234 g mol^-1
        Ionic concentration range: 0.09953429399586157 mol L^-1 - 5.30461986641539 mol L^-1
        Density range: 1.00116 g mL^-1 - 1.19412 g mL^-1
    ==================================================================
    """

    # Invalid cation/anion charges
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", cation_charge = -1,
    )
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", anion_charge = 1,
    )
    # Custom ion PDB file without an explicit charge
    @test_throws ArgumentError SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb", cation_pdbfile = "$test_dir/data/SOD.pdb",
    )

    # Neutral solute, no added salt: only water fills the box
    system = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 0)
    tmp_input_file = tempname()*".inp"
    rm(tmp_input_file, force=true)
    r = write_packmol_input(system; ionic_concentration = 0.0, margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r[1] > 0  # water
    @test r[2] == 0  # cations
    @test r[3] == 0  # anions
    @test r[4] ≈ [118.81, 118.81, 118.81]u"Å"
    input_text = read(tmp_input_file, String)
    @test occursin("structure $(Packmol.RecipesDirectory)/test/data/poly_h.pdb", input_text)
    @test !occursin("structure $(system.cation_pdbfile)", input_text)
    @test !occursin("structure $(system.anion_pdbfile)", input_text)

    # Neutral solute, physiological ionic strength: equal numbers of cations and anions
    rm(tmp_input_file, force=true)
    r = write_packmol_input(system; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test isfile(tmp_input_file)
    @test r[2] == r[3]
    @test r[2] > 0
    input_text = read(tmp_input_file, String)
    @test occursin("structure $(system.cation_pdbfile)", input_text)
    @test occursin("structure $(system.anion_pdbfile)", input_text)

    # The generated input file must be valid Packmol syntax for the native engine
    # itself, not just the (more lenient) legacy Fortran binary — regression test
    # for a stray comma once present in the `pbc` line, which the native parser
    # rejected outright.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == 1 + r[1] + r[2] + r[3]
    rm(tmp_input_file, force=true)

    # Charged solute (+5): extra anions neutralize it, on top of the bulk salt
    system_charged = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = 5)
    rm(tmp_input_file, force=true)
    r_charged = write_packmol_input(system_charged; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test r_charged[3] - r_charged[2] == 5  # 5 extra anions relative to cations

    # Negatively charged solute (-3): extra cations neutralize it
    system_neg = SolutionBoxUWI(solute_pdbfile = "$test_dir/data/poly_h.pdb", solute_charge = -3)
    rm(tmp_input_file, force=true)
    r_neg = write_packmol_input(system_neg; ionic_concentration = 0.15u"mol/L", margin = 20.0, cubic = true, input = tmp_input_file, debug = true)
    @test r_neg[2] - r_neg[3] == 3  # 3 extra cations relative to anions

    # Custom divalent cation (charge +2) neutralizing a negative solute charge: an odd
    # magnitude cannot be exactly cancelled by ions of charge magnitude 2
    system_divalent = SolutionBoxUWI(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solute_charge = -3,
        cation_pdbfile = "$test_dir/data/SOD.pdb",
        cation_charge = 2,
    )
    @test_throws ArgumentError write_packmol_input(
        system_divalent; ionic_concentration = 0.15u"mol/L", margin = 20.0, input = tmp_input_file,
    )

    rm(tmp_input_file, force=true)
end
