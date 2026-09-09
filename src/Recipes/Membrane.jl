#
# Membrane: a lipid bilayer or monolayer, solvated above (and, for a
# bilayer, also below) the membrane, under orthorhombic periodic boundary
# conditions.
#
# The membrane normal is always the z axis. Each lipid type is characterized
# by a single "head" atom and a single "tail" atom (1-based indices into its
# own PDB file); `h_i = norm(head_i - tail_i)`, the Euclidean distance
# between them in the given (arbitrarily oriented) PDB template, is taken as
# that lipid's intrinsic head-to-tail length — this is what the packing
# constraints below try to stretch each copy of that lipid to, regardless of
# how it happens to be rotated at any point during optimization.
#
export Membrane

# One lipid molecule's placement within one leaflet: `zlo`/`zhi` bound the
# whole molecule (all atoms), and the head/tail atoms are additionally
# pinned within `flex` Å of whichever edge (`head_low` selects which) faces
# this leaflet's solvent-facing anchor plane. All z values are already
# shifted so that z=0 is the box center (matching `_recipe_unitcell`'s
# `unitcell_center = (0,0,0)`).
struct _LipidPlacement
    pdbfile::String
    number::Int
    head_idx::Int
    tail_idx::Int
    zlo::Float64
    zhi::Float64
    head_low::Bool
    flex::Float64
end

# One solvent molecule's placement within one slab: `zlo`/`zhi` bound the
# whole molecule. Also already shifted to be centered at z=0.
struct _SolventPlacement
    pdbfile::String
    number::Int
    zlo::Float64
    zhi::Float64
end

mutable struct Membrane <: Recipe
    type::Symbol
    lipid_pdbfiles::Vector{String}
    lipid_head::Vector{Int}
    lipid_tail::Vector{Int}
    lipid_head_tail_length::Vector{typeof(1.0u"Å")}
    lipid_weight::Vector{Float64}
    lipid_molar_mass::Vector{typeof(1.0u"g/mol")}
    area_per_lipid::typeof(1.0u"Å^2")
    total_area::Union{Nothing,typeof(1.0u"Å^2")}
    total_lipids::Union{Nothing,Int}
    flexibility::Float64
    solvent_pdbfiles::Vector{String}
    solvent_weight::Vector{Float64}
    solvent_molar_mass::Vector{typeof(1.0u"g/mol")}
    solvent_layer_width::typeof(1.0u"Å")
    solvent_density::typeof(1.0u"g/mL")
end

"""
    Membrane(;
        type::Symbol = :bilayer, # or :monolayer
        lipids::Union{String,Vector{String}},
        lipid_head::Vector{<:Integer},
        lipid_tail::Vector{<:Integer},
        lipid_molar_ratio::Vector{<:Real},
        lipid_molar_mass = nothing, # optional
        area_per_lipid::Number,
        total_area::Union{Nothing,Number} = nothing, # xor total_lipids
        total_lipids::Union{Nothing,Integer} = nothing, # xor total_area
        flexibility::Real = 0.25,
        solvent::Union{String,Vector{String}},
        solvent_molar_ratio = nothing, # optional, defaults to equal weights
        solvent_molar_mass = nothing, # optional
        solvent_layer_width::Number,
        solvent_density::Number,
    )

Setup a lipid bilayer or monolayer (`type`), optionally with more than one lipid
(mixed at `lipid_molar_ratio`), solvated by one or more solvents (mixed at
`solvent_molar_ratio`). The system is always built with orthorhombic periodic
boundary conditions, with the membrane normal along z.

For each lipid, `lipid_head`/`lipid_tail` are the 1-based indices, into that lipid's
own PDB file, of a single atom marking the polar head and a single atom marking the
end of the tail. The Euclidean distance between them (in the given, arbitrarily
oriented, PDB template) is taken as that lipid's intrinsic head-to-tail length; the
membrane's own thickness is set by the *longest* such length among all the lipids
given (so a shorter lipid mixed in with longer ones simply doesn't reach all the way
to the bilayer's midplane, which is physically reasonable). Each lipid is packed
with two whole-molecule `plane` constraints bounding it to its leaflet's slab, plus
two per-atom `plane` constraints (applied via `atoms <idx> ... end atoms`, exactly as
in a hand-written Packmol input file) pinning its head atom near the leaflet's
solvent-facing edge and its tail atom near the leaflet's inner edge, each within
`flexibility` (a fraction, `0 < flexibility < 0.5`, of that lipid's own head-to-tail
length) of the exact edge — this is what gives the packing "room" to tilt/wobble each
lipid instead of forcing it rigidly upright.

For `type = :bilayer`, solvent fills the box above the top leaflet's heads and below
the bottom leaflet's heads (both leaflets share the same composition/counts); for
`type = :monolayer`, solvent fills the box only above the single leaflet's heads,
and the tails point towards the opposite (bottom) box face.

Exactly one of `area_per_lipid` or `total_area` and exactly one of ... (see below)
must be given: `area_per_lipid` (Å², the average lateral area per lipid molecule)
is always required; together with exactly one of `total_area` (Å², the box's lateral
footprint) or `total_lipids` (the total lipid count, summed over both leaflets for a
bilayer) it fixes both the box's lateral size and the lipid count per leaflet — the
one not given is computed from the other, mirroring how `box_sides`/`margin` work in
`SolutionBoxUS`.

`solvent_layer_width` (Å) is the thickness of each solvent slab, and `solvent_density`
(`g/mL` by default) is the density of the solvent mixture, used together with
`solvent_molar_ratio` to size the number of molecules of each solvent species (as in
`SolutionBoxUSC`'s cossolvent mixture, but with a directly given density instead of
a density table).

Molar masses are computed from the atom types in the PDB files if not provided.

"""
function Membrane(;
    type::Symbol=:bilayer,
    lipids::Union{String,Vector{String}},
    lipid_head::Vector{<:Integer},
    lipid_tail::Vector{<:Integer},
    lipid_molar_ratio::Vector{<:Real},
    lipid_molar_mass::Union{Nothing,Vector{<:Number}}=nothing,
    area_per_lipid::Number,
    total_area::Union{Nothing,Number}=nothing,
    total_lipids::Union{Nothing,Integer}=nothing,
    flexibility::Real=0.25,
    solvent::Union{String,Vector{String}},
    solvent_molar_ratio::Union{Nothing,Vector{<:Real}}=nothing,
    solvent_molar_mass::Union{Nothing,Vector{<:Number}}=nothing,
    solvent_layer_width::Number,
    solvent_density::Number,
)
    type in (:bilayer, :monolayer) || throw(ArgumentError("type must be :bilayer or :monolayer, got :$type"))

    lipids isa String && (lipids = [lipids])
    nlip = length(lipids)
    (length(lipid_head) == nlip && length(lipid_tail) == nlip && length(lipid_molar_ratio) == nlip) ||
        throw(ArgumentError("lipids, lipid_head, lipid_tail, and lipid_molar_ratio must all have the same length."))
    !isnothing(lipid_molar_mass) && length(lipid_molar_mass) != nlip &&
        throw(ArgumentError("lipid_molar_mass must have the same length as lipids."))
    all(>(0), lipid_molar_ratio) || throw(ArgumentError("lipid_molar_ratio must be all positive."))

    isnothing(total_area) == isnothing(total_lipids) &&
        throw(ArgumentError("Exactly one of total_area or total_lipids must be provided."))
    !isnothing(total_lipids) && total_lipids <= 0 &&
        throw(ArgumentError("total_lipids must be positive, got $total_lipids."))

    0 < flexibility < 0.5 || throw(ArgumentError("flexibility must satisfy 0 < flexibility < 0.5, got $flexibility."))

    lipid_atoms = [read_pdb(f) for f in lipids]
    for i in 1:nlip
        natoms_i = length(lipid_atoms[i])
        (1 <= lipid_head[i] <= natoms_i) ||
            throw(ArgumentError("lipid_head[$i] = $(lipid_head[i]) is out of range for $(lipids[i]) ($natoms_i atoms)."))
        (1 <= lipid_tail[i] <= natoms_i) ||
            throw(ArgumentError("lipid_tail[$i] = $(lipid_tail[i]) is out of range for $(lipids[i]) ($natoms_i atoms)."))
        lipid_head[i] == lipid_tail[i] &&
            throw(ArgumentError("lipid_head[$i] and lipid_tail[$i] must refer to different atoms."))
    end
    lipid_head_tail_length = [
        begin
            h, t = lipid_atoms[i][lipid_head[i]], lipid_atoms[i][lipid_tail[i]]
            sqrt((h.x - t.x)^2 + (h.y - t.y)^2 + (h.z - t.z)^2) * u"Å"
        end
        for i in 1:nlip
    ]
    lipid_molar_mass = isnothing(lipid_molar_mass) ? [mass(atoms) for atoms in lipid_atoms] : collect(lipid_molar_mass)
    lipid_molar_mass = _ensure_unit.(lipid_molar_mass, u"g/mol")
    lipid_weight = lipid_molar_ratio ./ sum(lipid_molar_ratio)

    area_per_lipid = _ensure_unit(area_per_lipid, u"Å^2")
    ustrip(area_per_lipid) > 0 || throw(ArgumentError("area_per_lipid must be positive."))
    !isnothing(total_area) && (total_area = _ensure_unit(total_area, u"Å^2"))

    solvent isa String && (solvent = [solvent])
    nsolv = length(solvent)
    solvent_molar_ratio = isnothing(solvent_molar_ratio) ? fill(1.0, nsolv) : solvent_molar_ratio
    length(solvent_molar_ratio) == nsolv ||
        throw(ArgumentError("solvent_molar_ratio must have the same length as solvent."))
    !isnothing(solvent_molar_mass) && length(solvent_molar_mass) != nsolv &&
        throw(ArgumentError("solvent_molar_mass must have the same length as solvent."))
    all(>(0), solvent_molar_ratio) || throw(ArgumentError("solvent_molar_ratio must be all positive."))
    solvent_molar_mass = isnothing(solvent_molar_mass) ?
        [mass(read_pdb(f)) for f in solvent] : collect(solvent_molar_mass)
    solvent_molar_mass = _ensure_unit.(solvent_molar_mass, u"g/mol")
    solvent_weight = solvent_molar_ratio ./ sum(solvent_molar_ratio)

    solvent_layer_width = _ensure_unit(solvent_layer_width, u"Å")
    ustrip(solvent_layer_width) > 0 || throw(ArgumentError("solvent_layer_width must be positive."))
    if unit(solvent_density) == NoUnits
        @warn "Density units not provided, assuming g/mL." _file=nothing _line=nothing
        solvent_density = solvent_density * 1.0u"g/mL"
    end
    ustrip(solvent_density) > 0 || throw(ArgumentError("solvent_density must be positive."))

    return Membrane(
        type, collect(lipids), collect(lipid_head), collect(lipid_tail),
        lipid_head_tail_length, lipid_weight, lipid_molar_mass,
        area_per_lipid, total_area, total_lipids, Float64(flexibility),
        collect(solvent), solvent_weight, solvent_molar_mass,
        solvent_layer_width, solvent_density,
    )
end

function Base.show(io::IO, ::MIME"text/plain", system::Membrane)
    print(io, chomp("""
    ==================================================================
    Membrane properties ($(system.type)):
    ==================================================================
        Lipids: $(join(basename.(system.lipid_pdbfiles), ", "))
        Lipid molar ratio (normalized): $(round.(system.lipid_weight; digits=4))
        Lipid head-to-tail lengths: $(join(system.lipid_head_tail_length, ", "))
        Solvents: $(join(basename.(system.solvent_pdbfiles), ", "))
        Solvent molar ratio (normalized): $(round.(system.solvent_weight; digits=4))
        Solvent density: $(system.solvent_density)
        Solvent layer width: $(system.solvent_layer_width)
        Area per lipid: $(system.area_per_lipid)
        $(isnothing(system.total_lipids) ? "Total area: $(system.total_area)" : "Total lipids: $(system.total_lipids)")
        Flexibility: $(system.flexibility)
    ==================================================================
    """))
end

#
# Shared box-sizing/molecule-count computation for Membrane, used by both
# `write_packmol_input` (which prints/writes `summary` to a `.inp` file) and
# `packmol` (which builds a `PackmolSystem` directly and never writes one).
#
function _setup(system::Membrane)
    (; type, lipid_pdbfiles, lipid_head, lipid_tail, lipid_head_tail_length, lipid_weight, lipid_molar_mass,
       area_per_lipid, total_area, total_lipids, flexibility,
       solvent_pdbfiles, solvent_weight, solvent_molar_mass, solvent_layer_width, solvent_density) = system

    n_leaflets = type == :bilayer ? 2 : 1
    nlip = length(lipid_pdbfiles)
    nsolv = length(solvent_pdbfiles)
    d = maximum(lipid_head_tail_length) # longest lipid: sets the membrane's nominal thickness

    if isnothing(total_lipids)
        lipids_per_leaflet = round(Int, ustrip(u"Å^2", total_area) / ustrip(u"Å^2", area_per_lipid))
        total_area_used = total_area
    else
        lipids_per_leaflet = round(Int, total_lipids / n_leaflets)
        total_area_used = lipids_per_leaflet * area_per_lipid
    end
    lipids_per_leaflet > 0 || throw(ArgumentError(
        "Computed zero lipids per leaflet: increase total_lipids/total_area, or decrease area_per_lipid."
    ))
    total_lipids_actual = lipids_per_leaflet * n_leaflets

    Lxy = sqrt(ustrip(u"Å^2", total_area_used)) * u"Å"

    # Per-lipid-type count for one leaflet (both leaflets share it, for a
    # bilayer), corrected so the counts sum exactly to `lipids_per_leaflet`.
    lipid_counts = round.(Int, lipid_weight .* lipids_per_leaflet)
    lipid_counts[end] += lipids_per_leaflet - sum(lipid_counts)

    # Membrane/box z geometry (Å, unshifted: z=0 at the very bottom of the box)
    w = ustrip(u"Å", solvent_layer_width)
    dÅ = ustrip(u"Å", d)
    Lz = type == :bilayer ? 2w + 2dÅ : w + dÅ

    # Leaflet anchor plane (the z where every lipid type's head sits, shared
    # across lipid types since they all face the same solvent interface) and
    # `head_low` (whether the head sits at the low-z or high-z edge of the
    # leaflet's own slab), one pair per leaflet.
    anchors, head_lows = type == :bilayer ? ([w, w + 2dÅ], (true, false)) : ([dÅ], (false,))

    lipid_placements = _LipidPlacement[]
    for (A, head_low) in zip(anchors, head_lows)
        for i in 1:nlip
            lipid_counts[i] == 0 && continue
            hi = ustrip(u"Å", lipid_head_tail_length[i])
            flex = flexibility * hi
            zlo, zhi = head_low ? (A, A + hi) : (A - hi, A)
            push!(lipid_placements, _LipidPlacement(
                lipid_pdbfiles[i], lipid_counts[i], lipid_head[i], lipid_tail[i],
                zlo - Lz / 2, zhi - Lz / 2, head_low, flex,
            ))
        end
    end

    # Solvent slab z ranges (one slab per leaflet: 2 for a bilayer, 1 for a
    # monolayer), all slabs sharing the same composition/counts.
    slab_bounds = type == :bilayer ? [(0.0, w), (w + 2dÅ, Lz)] : [(dÅ, Lz)]
    mean_M = sum(solvent_weight[k] * solvent_molar_mass[k] for k in 1:nsolv)
    slab_volume = Lxy^2 * solvent_layer_width
    total_mass = uconvert(u"g", solvent_density * slab_volume)
    total_moles = total_mass / mean_M
    solvent_counts = [round(Int, ustrip(solvent_weight[k] * total_moles * Unitful.Na)) for k in 1:nsolv]

    solvent_placements = _SolventPlacement[]
    for (zlo_raw, zhi_raw) in slab_bounds
        for k in 1:nsolv
            solvent_counts[k] == 0 && continue
            push!(solvent_placements, _SolventPlacement(
                solvent_pdbfiles[k], solvent_counts[k], zlo_raw - Lz / 2, zhi_raw - Lz / 2,
            ))
        end
    end

    unitcell = Matrix{Float64}(Diagonal([ustrip(u"Å", Lxy), ustrip(u"Å", Lxy), Lz]))

    lipid_lines = join(
        ("    $(basename(lipid_pdbfiles[i])): $(lipid_counts[i]) per leaflet, " *
         "head-to-tail length = $(lipid_head_tail_length[i])" for i in 1:nlip),
        "\n",
    )
    solvent_lines = join(
        ("    $(basename(solvent_pdbfiles[k])): $(solvent_counts[k]) per slab" for k in 1:nsolv),
        "\n",
    )
    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Membrane type = $type ($n_leaflets leaflet(s))
        Membrane (nominal) thickness = $(n_leaflets * d)
        Box lateral area = $total_area_used
        Periodic box = $(_unitcell_description(unitcell))

        Lipids (per leaflet, $lipids_per_leaflet total, $total_lipids_actual overall):
        $lipid_lines

        Solvents (per slab):
        $solvent_lines

        ==================================================================
        """
    return (; lipid_placements, solvent_placements, unitcell, summary)
end

"""
    write_packmol_input(
        system::Membrane;
        input="membrane.inp",
        output="membrane.pdb",
    )

Function that generates an input file for Packmol for a lipid membrane system.

"""
function write_packmol_input(
    system::Membrane;
    input="membrane.inp",
    output="membrane.pdb",
    # testing option
    debug=false,
)
    (; lipid_placements, solvent_placements, unitcell, summary) = _setup(system)
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
            """)
        for lp in lipid_placements
            head_line = lp.head_low ?
                "below plane 0. 0. 1. $(lp.zlo + lp.flex)" : "above plane 0. 0. 1. $(lp.zhi - lp.flex)"
            tail_line = lp.head_low ?
                "above plane 0. 0. 1. $(lp.zhi - lp.flex)" : "below plane 0. 0. 1. $(lp.zlo + lp.flex)"
            println(io,
                """
                structure $(lp.pdbfile)
                    number $(lp.number)
                    above plane 0. 0. 1. $(lp.zlo)
                    below plane 0. 0. 1. $(lp.zhi)
                    atoms $(lp.head_idx)
                        $head_line
                    end atoms
                    atoms $(lp.tail_idx)
                        $tail_line
                    end atoms
                end structure
                """)
        end
        for sp in solvent_placements
            println(io,
                """
                structure $(sp.pdbfile)
                    number $(sp.number)
                    above plane 0. 0. 1. $(sp.zlo)
                    below plane 0. 0. 1. $(sp.zhi)
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
        return lipid_placements, solvent_placements, unitcell
    else
        return nothing
    end
end # function write_packmol_input

"""
    packmol(
        system::Membrane;
        output="membrane.pdb",
        kwargs...,
    )

Builds and packs a lipid membrane system directly, entirely in memory: equivalent to
calling [`write_packmol_input`](@ref write_packmol_input(::Membrane)) followed by
`packmol` on the resulting file, except no `.inp` file is ever written.

`output` behaves as in `write_packmol_input`. Any other keyword (`nloop`, `iprint`,
`seed`, `optimizer`, ...) is forwarded to the packing engine — see
`packmol(::PackmolSystem)`.

Returns the built `PackmolSystem`, with the packing outcome in its `.status` field —
see `packmol(::PackmolSystem)`.

"""
function packmol(
    system::Membrane;
    output="membrane.pdb",
    kwargs...,
)
    (; lipid_placements, solvent_placements, unitcell) = _setup(system)
    structure_types = StructureType{3,Float64}[]
    for lp in lipid_placements
        atoms = read_pdb(lp.pdbfile)
        natoms = length(atoms)
        reference_coordinates = [SVector{3,Float64}(a.x, a.y, a.z) for a in atoms]
        radii = fill(1.0, natoms)
        c1 = AbovePlane([0.0, 0.0, 1.0], lp.zlo)
        c2 = BelowPlane([0.0, 0.0, 1.0], lp.zhi)
        head_c = lp.head_low ?
            BelowPlane([0.0, 0.0, 1.0], lp.zlo + lp.flex) : AbovePlane([0.0, 0.0, 1.0], lp.zhi - lp.flex)
        tail_c = lp.head_low ?
            AbovePlane([0.0, 0.0, 1.0], lp.zhi - lp.flex) : BelowPlane([0.0, 0.0, 1.0], lp.zlo + lp.flex)
        atom_constraints = [Int[1, 2] for _ in 1:natoms]
        push!(atom_constraints[lp.head_idx], 3)
        push!(atom_constraints[lp.tail_idx], 4)
        push!(structure_types, StructureType{3,Float64}(;
            filename=lp.pdbfile, natoms, atoms, number_of_molecules=lp.number,
            reference_coordinates, radii,
            constraints=AnyConstraint{Float64}[c1, c2, head_c, tail_c],
            atom_constraints,
        ))
    end
    for sp in solvent_placements
        push!(structure_types, structure_type(sp.pdbfile; number=sp.number,
            constraints=[AbovePlane([0.0, 0.0, 1.0], sp.zlo), BelowPlane([0.0, 0.0, 1.0], sp.zhi)],
        ))
    end
    packmol_system = PackmolSystem(structure_types;
        output, tolerance=2.0, add_box_sides=true, seed=-1, _recipe_unitcell(unitcell)...,
    )
    return packmol(packmol_system; kwargs...)
end
