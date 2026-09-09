mutable struct SolutionBoxUS <: Recipe
    solute_pdbfile::String
    solvent_pdbfile::String
    density::Quantity
    solute_molar_mass::Quantity
    solvent_molar_mass::Quantity
end

#
# Builds the periodic unit cell (a plain 3×3 `Matrix{Float64}`, in Å) enclosing
# the solute plus the requested box sides/margin, in the shape requested by
# `pbc` (`:cubic`, `:orthorhombic`, `:dodecahedral`, or `:octahedral`):
#   - `:orthorhombic` uses `box_sides` (or `solute_extrema .+ 2margin`) as given,
#     possibly with unequal sides.
#   - `:cubic` forces all three sides to their maximum, like `:orthorhombic`
#     but with a cube instead of a general box.
#   - `:dodecahedral`/`:octahedral` build a rhombic dodecahedron/truncated
#     octahedron (see `dodecahedral_unitcell`/`octahedral_unitcell`) of size
#     `d = maximum(box_sides)` — the same "take the largest requested side"
#     logic as `:cubic`, just for that shape instead of a cube.
#
function set_unitcell(system, box_sides, margin, pbc::Symbol)
    pbc in (:cubic, :orthorhombic, :dodecahedral, :octahedral) || throw(ArgumentError(
        "pbc must be :cubic, :orthorhombic, :dodecahedral, or :octahedral, got :$pbc"
    ))
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
    unitcell = if pbc == :dodecahedral
        d = ustrip(u"Å", maximum(box_sides))
        dodecahedral_unitcell(Float64, d)
    elseif pbc == :octahedral
        d = ustrip(u"Å", maximum(box_sides))
        octahedral_unitcell(Float64, d)
    else
        pbc == :cubic && (box_sides = fill(maximum(box_sides), 3))
        Matrix{Float64}(Diagonal(ustrip.(u"Å", box_sides)))
    end
    return unitcell, solute_extrema
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

#
# Shared box-sizing/molecule-count computation for SolutionBoxUS, used by both
# `write_packmol_input` (which prints/writes `summary` to a `.inp` file) and
# `packmol` (which builds a `PackmolSystem` directly and never writes one).
#
function _setup(
    system::SolutionBoxUS,
    box_sides::Union{AbstractVector{<:Number},Nothing},
    margin::Union{<:Number,Nothing},
    pbc::Symbol,
)
    (; solvent_pdbfile, solute_molar_mass, solvent_molar_mass) = system

    # molar masses (g/mol)
    Mp = solute_molar_mass
    Mw = solvent_molar_mass

    # Density of pure solvent (g/mL)
    ρs = system.density

    # Molarity of pure solvent (mol/L)
    ms = uconvert(u"mol/L", ρs / Mw)

    # Convert solvent concentration in molecules/Å³
    cs = cconvert(ms, "mol/L" => "molecules/Å^-3")

    # Set unit cell and volume
    unitcell, solute_extrema = set_unitcell(system, box_sides, margin, pbc)
    vbox = det(unitcell) * u"Å^3"

    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the pure solvent
    vs = vbox - uconvert(u"Å^3", Mp / ρs / Unitful.Na)

    # number of solvent molecules (molecules/Å³ * Å³)
    ns = round(Int, cs * vs)

    # Number of solvent molecules
    ns == 0 && throw(ArgumentError("Number of solvent molecules is zero."))

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
        Periodic box (pbc = :$pbc) = $(_unitcell_description(unitcell))

        Solute molar mass = $Mp
        Solvent molar mass = $Mw

        Number of solvent ($(basename(solvent_pdbfile))) molecules = $ns

        ==================================================================
        """
    return (; ns, unitcell, summary)
end

"""
    write_packmol_input(
        system::SolutionBoxUS;
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        pbc::Symbol = :cubic,
    )

Function that generates an input file for Packmol for a Solute + Solvent system.

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternatively, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

`pbc` selects the shape of the periodic cell: `:cubic` (the default) forces all 3
sides to their maximum, `:orthorhombic` keeps `box_sides`/`margin` as given (possibly
with unequal sides), and `:dodecahedral`/`:octahedral` build a rhombic
dodecahedron/truncated octahedron cell (see [`dodecahedral_unitcell`](@ref)/
[`octahedral_unitcell`](@ref)) of size equal to that same maximum side.

"""
function write_packmol_input(
    system::SolutionBoxUS;
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing} = nothing,
    margin::Union{<:Number,Nothing} = nothing,
    pbc::Symbol = :cubic,
    # testing option
    debug = false,
)
    (; solute_pdbfile, solvent_pdbfile) = system
    (; ns, unitcell, summary) = _setup(system, box_sides, margin, pbc)
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
        a, b, c, α, β, γ = _unitcell_abc_angles(unitcell)
        println(io,
            """
            #
            tolerance 2.0
            output $output
            add_box_sides 1.0
            filetype pdb
            seed -1
            packall
            unitcell $a $b $c $α $β $γ

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
        a, b, c, = _unitcell_abc_angles(unitcell)
        return ns, [a, b, c] * u"Å"
    else
        return nothing
    end
end # function write_packmol_input

"""
    packmol(
        system::SolutionBoxUS;
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Number}, # or
        margin::Number,
        pbc::Symbol = :cubic,
        kwargs...,
    )

Builds and packs a Solute + Solvent system directly, entirely in memory: equivalent to
calling [`write_packmol_input`](@ref write_packmol_input(::SolutionBoxUS)) followed by
`packmol` on the resulting file, except no `.inp` file is ever written.

`output`, `box_sides`, `margin`, and `pbc` behave as in `write_packmol_input`. Any other
keyword (`nloop`, `iprint`, `seed`, `optimizer`, ...) is forwarded to the packing engine —
see `packmol(::PackmolSystem)`.

Returns the built `PackmolSystem`, with the packing outcome in its `.status` field —
see `packmol(::PackmolSystem)`.

"""
function packmol(
    system::SolutionBoxUS;
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Number},Nothing}=nothing,
    margin::Union{<:Number,Nothing}=nothing,
    pbc::Symbol=:cubic,
    kwargs...,
)
    (; solute_pdbfile, solvent_pdbfile) = system
    (; ns, unitcell) = _setup(system, box_sides, margin, pbc)
    structure_types = [
        _fixed_solute_structure_type(solute_pdbfile),
        structure_type(solvent_pdbfile; number=ns),
    ]
    packmol_system = PackmolSystem(structure_types;
        output, tolerance=2.0, add_box_sides=true, seed=-1, _recipe_unitcell(unitcell)...,
    )
    return packmol(packmol_system; kwargs...)
end
