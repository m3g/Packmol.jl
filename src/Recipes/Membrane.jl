#
# Membrane: a lipid bilayer or monolayer, solvated above (and, for a bilayer, also
# below) the membrane, under orthorhombic periodic boundary conditions, with the
# membrane normal along z.
#
# Each lipid type is characterized by a head atom and its chain-end atoms (1-based
# indices into its own PDB file); its length L is the distance from the head atom to
# the midpoint of the chain ends. The chain ends are kept within `_CHAIN_END_PLANE` of
# the middle of the membrane, and the head above `_CHAIN_END_PLANE + (1 - flexibility) L`
# from it: a narrow band of the heads keeps the lipids upright, which is what makes
# them pack (with `flexibility = 0.5`, most packings of the lipids of the `lipids/`
# directory fail). No constraint bounds the whole molecule.
#
export Membrane

# The built-in lipids (`lipids/`, see its README): head atom, chain-end atoms, and the
# area per lipid of the actual membrane (Å²)
const _MEMBRANE_LIPIDS_DIR = joinpath(@__DIR__, "lipids")
const _BUILTIN_LIPIDS = Dict(
    "POPC" => (head=1, chain_ends=[88, 131], area=64.0),
    "POPE" => (head=1, chain_ends=[79, 122], area=57.0),
    "POPS" => (head=1, chain_ends=[81, 124], area=63.0),
    "POPG" => (head=1, chain_ends=[81, 124], area=66.0),
    "DPPC" => (head=1, chain_ends=[84, 127], area=63.0),
    "DOPC" => (head=1, chain_ends=[88, 135], area=67.0),
    "DMPC" => (head=1, chain_ends=[78, 115], area=61.0),
)
# The default area per lipid is this much larger than that of the actual membrane: at the
# actual areas, a packing of the rigid conformations gets stuck in some runs (DMPC, POPE,
# and POPC/POPS 4:1 in both of two runs), at 5% more all of them pack in seconds.
const _MEMBRANE_AREA_FACTOR = 1.05
# Distance (Å) of the plane of the chain ends from the middle of a bilayer (from the
# bottom of a monolayer), and of the solvent from the plane of the heads
const _CHAIN_END_PLANE = 6.0
const _SOLVENT_GAP = 2.0

# One lipid type in one leaflet: its head atom above `head_plane` and its chain ends
# below `tail_plane` (`upper`), or the mirror (lower leaflet of a bilayer); `floor`, for
# a monolayer, bounds the whole molecule from below. z = 0 is the center of the box.
struct _LipidPlacement
    pdbfile::String
    number::Int
    head_idx::Int
    tail_idx::Vector{Int}
    head_plane::Float64
    tail_plane::Float64
    upper::Bool
    floor::Union{Nothing,Float64}
end

# One solvent in one slab: `zlo`/`zhi` bound the whole molecule. z = 0 is the center of
# the box.
struct _SolventPlacement
    pdbfile::String
    number::Int
    zlo::Float64
    zhi::Float64
end

mutable struct Membrane <: Recipe
    type::Symbol
    lipid_names::Vector{String}
    lipid_pdbfiles::Vector{String}
    lipid_head::Vector{Int}
    lipid_tail::Vector{Vector{Int}}
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

_mandatory_keywords(::Type{Membrane}) = (:lipids,)

# The built-in lipid of `name` (case insensitive), or `nothing`
_builtin_lipid(name::AbstractString) = get(_BUILTIN_LIPIDS, uppercase(name), nothing)
_is_builtin_lipid(name::AbstractString) = !isfile(expanduser(name)) && !isnothing(_builtin_lipid(name))

"""
    Membrane(;
        lipids::Union{String,Vector{String}},
        lipid_molar_ratio = nothing, # equal amounts if not given
        total_lipids::Union{Nothing,Integer} = nothing, # xor total_area
        total_area::Union{Nothing,Number} = nothing, # Å², xor total_lipids
        type::Symbol = :bilayer, # or :monolayer
        area_per_lipid = nothing, # Å², automatic for the built-in lipids
        lipid_head = nothing, # automatic for the built-in lipids
        lipid_tail = nothing, # automatic for the built-in lipids
        lipid_molar_mass = nothing,
        flexibility::Real = 0.15,
        solvent = nothing, # water
        solvent_molar_ratio = nothing, # equal amounts if not given
        solvent_molar_mass = nothing,
        solvent_layer_width::Number = 15.0, # Å
        solvent_density::Number = 1.0u"g/mL",
    )

Setup a lipid bilayer or monolayer (`type`) of one or more lipids, mixed at
`lipid_molar_ratio`, solvated by water or by a mixture of solvents. The system is
built with orthorhombic periodic boundary conditions, with the membrane normal along z.

The lipids are the names of built-in lipids, `"POPC"`, `"POPE"`, `"POPS"`, `"POPG"`,
`"DPPC"`, `"DOPC"` and `"DMPC"`, or PDB files. For the built-in lipids, only the
composition and the size of the membrane are needed: for example,

```julia
Membrane(lipids = ["POPC", "POPE"], lipid_molar_ratio = [3, 1], total_lipids = 200)
```

is a POPC/POPE 3:1 bilayer of 200 lipids in water.

The size of the membrane is given by exactly one of `total_lipids` (the number of lipids,
summed over both leaflets of a bilayer) or `total_area` (Å², the lateral area of the
box), with `area_per_lipid` (Å²): by default (built-in lipids only), 5% more than the
molar average of the areas per lipid of the actual membranes (POPC 64, POPE 57, POPS 63,
POPG 66, DPPC 63, DOPC 67, DMPC 61 Å²). At the actual areas, the packing of rigid
lipids gets stuck in some runs; a molecular dynamics simulation of the system then
equilibrates the area. Both leaflets of a bilayer have the same composition.

For a lipid given as a PDB file, `lipid_head` (the index of an atom of the polar head in
the file) and `lipid_tail` (the indices of the atoms at the ends of its chains, a vector
for each lipid, or a single index) must be given for all the lipids, in the order of
`lipids`, and `area_per_lipid` too. Packing is fastest with straight, compact
conformations, with the chains side by side, as the built-in ones.

The length L of each lipid is the distance from its head atom to the midpoint of its
chain ends. The chain ends are packed within 6 Å of the middle of the bilayer (of the
bottom of a monolayer), and the head above 6 Å + (1 - `flexibility`) L: `flexibility`
(`0 < flexibility < 1`) sets how far from upright the lipids may be. Wider bands of the
heads (larger `flexibility`) make the packing slower, or make it fail.

The solvent fills slabs of `solvent_layer_width` (Å) beyond the heads (above and below a
bilayer, above a monolayer), of density `solvent_density` (`g/mL` by default), mixed at
`solvent_molar_ratio`. By default, it is water at 1 g/mL. Molar masses are computed from
the atom types in the PDB files if not given.

"""
function Membrane(;
    lipids::Union{String,Vector{String}} = _missing_keyword(Membrane, :lipids),
    lipid_molar_ratio::Union{Nothing,Vector{<:Real}}=nothing,
    total_lipids::Union{Nothing,Integer}=nothing,
    total_area::Union{Nothing,Number}=nothing,
    type::Symbol=:bilayer,
    area_per_lipid::Union{Nothing,Number}=nothing,
    lipid_head::Union{Nothing,Vector{<:Integer}}=nothing,
    lipid_tail::Union{Nothing,Vector}=nothing,
    lipid_molar_mass::Union{Nothing,Vector{<:Number}}=nothing,
    flexibility::Real=0.15,
    solvent::Union{Nothing,String,Vector{String}}=nothing,
    solvent_molar_ratio::Union{Nothing,Vector{<:Real}}=nothing,
    solvent_molar_mass::Union{Nothing,Vector{<:Number}}=nothing,
    solvent_layer_width::Number=15.0,
    solvent_density::Number=1.0u"g/mL",
)
    type in (:bilayer, :monolayer) || throw(ArgumentError("type must be :bilayer or :monolayer, got :$type"))

    lipids isa String && (lipids = [lipids])
    nlip = length(lipids)
    nlip > 0 || throw(ArgumentError("At least one lipid must be given."))
    builtin = _is_builtin_lipid.(lipids)
    lipid_names = [builtin[i] ? uppercase(lipids[i]) : basename(lipids[i]) for i in 1:nlip]
    lipid_pdbfiles = [
        builtin[i] ? joinpath(_MEMBRANE_LIPIDS_DIR, lowercase(lipids[i]) * ".pdb") : expanduser(lipids[i])
        for i in 1:nlip
    ]
    for i in 1:nlip
        builtin[i] || isfile(lipid_pdbfiles[i]) || throw(ArgumentError("""\n
            Lipid `$(lipids[i])` is neither a built-in lipid ($(join(sort!(collect(keys(_BUILTIN_LIPIDS))), ", ")))
            nor an existing PDB file.

        """))
    end

    # Head and chain-end atoms: given for all the lipids, or automatic (built-in lipids)
    if isnothing(lipid_head) != isnothing(lipid_tail)
        throw(ArgumentError("lipid_head and lipid_tail must be given together."))
    end
    if isnothing(lipid_head)
        all(builtin) || throw(ArgumentError("""\n
            lipid_head and lipid_tail must be given for lipids given as PDB files
            ($(join(lipids[.!builtin], ", "))).

        """))
        lipid_head = [_builtin_lipid(l).head for l in lipids]
        lipid_tail = [copy(_builtin_lipid(l).chain_ends) for l in lipids]
    else
        lipid_tail = [t isa Integer ? [Int(t)] : collect(Int, t) for t in lipid_tail]
    end
    (length(lipid_head) == nlip && length(lipid_tail) == nlip) ||
        throw(ArgumentError("lipid_head and lipid_tail must have one element for each lipid."))

    lipid_molar_ratio = isnothing(lipid_molar_ratio) ? fill(1.0, nlip) : lipid_molar_ratio
    length(lipid_molar_ratio) == nlip ||
        throw(ArgumentError("lipid_molar_ratio must have one element for each lipid."))
    all(>(0), lipid_molar_ratio) || throw(ArgumentError("lipid_molar_ratio must be all positive."))
    !isnothing(lipid_molar_mass) && length(lipid_molar_mass) != nlip &&
        throw(ArgumentError("lipid_molar_mass must have one element for each lipid."))

    isnothing(total_area) == isnothing(total_lipids) &&
        throw(ArgumentError("Exactly one of total_area or total_lipids must be provided."))
    !isnothing(total_lipids) && total_lipids <= 0 &&
        throw(ArgumentError("total_lipids must be positive, got $total_lipids."))

    0 < flexibility < 1 || throw(ArgumentError("flexibility must satisfy 0 < flexibility < 1, got $flexibility."))

    lipid_atoms = [read_pdb(f) for f in lipid_pdbfiles]
    for i in 1:nlip
        natoms_i = length(lipid_atoms[i])
        (1 <= lipid_head[i] <= natoms_i) ||
            throw(ArgumentError("lipid_head[$i] = $(lipid_head[i]) is out of range for $(lipids[i]) ($natoms_i atoms)."))
        isempty(lipid_tail[i]) && throw(ArgumentError("lipid_tail[$i] is empty."))
        all(t -> 1 <= t <= natoms_i, lipid_tail[i]) ||
            throw(ArgumentError("lipid_tail[$i] = $(lipid_tail[i]) is out of range for $(lipids[i]) ($natoms_i atoms)."))
        lipid_head[i] in lipid_tail[i] &&
            throw(ArgumentError("lipid_head[$i] and lipid_tail[$i] must refer to different atoms."))
    end
    lipid_head_tail_length = [
        begin
            h = lipid_atoms[i][lipid_head[i]]
            ends = lipid_atoms[i][lipid_tail[i]]
            n = length(ends)
            mx, my, mz = sum(a.x for a in ends) / n, sum(a.y for a in ends) / n, sum(a.z for a in ends) / n
            sqrt((h.x - mx)^2 + (h.y - my)^2 + (h.z - mz)^2) * u"Å"
        end
        for i in 1:nlip
    ]
    lipid_molar_mass = isnothing(lipid_molar_mass) ? [mass(atoms) for atoms in lipid_atoms] : collect(lipid_molar_mass)
    lipid_molar_mass = _ensure_unit.(lipid_molar_mass, u"g/mol")
    lipid_weight = lipid_molar_ratio ./ sum(lipid_molar_ratio)

    if isnothing(area_per_lipid)
        all(builtin) || throw(ArgumentError("""\n
            area_per_lipid must be given for lipids given as PDB files
            ($(join(lipids[.!builtin], ", "))).

        """))
        area_per_lipid = _MEMBRANE_AREA_FACTOR * sum(lipid_weight[i] * _builtin_lipid(lipids[i]).area for i in 1:nlip)
    end
    area_per_lipid = _ensure_unit(area_per_lipid, u"Å^2")
    ustrip(area_per_lipid) > 0 || throw(ArgumentError("area_per_lipid must be positive."))
    !isnothing(total_area) && (total_area = _ensure_unit(total_area, u"Å^2"))

    solvent = isnothing(solvent) ? [joinpath(_MEMBRANE_LIPIDS_DIR, "water.pdb")] :
        solvent isa String ? [expanduser(solvent)] : expanduser.(solvent)
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
        type, lipid_names, lipid_pdbfiles, collect(Int, lipid_head), lipid_tail,
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
        Lipids: $(join(system.lipid_names, ", "))
        Lipid molar ratio (normalized): $(round.(system.lipid_weight; digits=4))
        Lipid head-to-tail lengths: $(join(round.(u"Å", system.lipid_head_tail_length; digits=2), ", "))
        Solvents: $(join(basename.(system.solvent_pdbfiles), ", "))
        Solvent molar ratio (normalized): $(round.(system.solvent_weight; digits=4))
        Solvent density: $(system.solvent_density)
        Solvent layer width: $(system.solvent_layer_width)
        Area per lipid: $(round(u"Å^2", system.area_per_lipid; digits=2))
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
    (; type, lipid_names, lipid_pdbfiles, lipid_head, lipid_tail, lipid_head_tail_length, lipid_weight,
       area_per_lipid, total_area, total_lipids, flexibility,
       solvent_pdbfiles, solvent_weight, solvent_molar_mass, solvent_layer_width, solvent_density) = system

    n_leaflets = type == :bilayer ? 2 : 1
    nlip = length(lipid_pdbfiles)
    nsolv = length(solvent_pdbfiles)

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

    # Planes, relative to the middle of the bilayer (the bottom of the monolayer): the
    # chain ends below _CHAIN_END_PLANE, the head of each lipid above its own plane, and
    # the solvent from _SOLVENT_GAP above the highest of them
    L = ustrip.(u"Å", lipid_head_tail_length)
    head_planes = [_CHAIN_END_PLANE + (1 - flexibility) * L[i] for i in 1:nlip]
    zsolvent = maximum(head_planes) + _SOLVENT_GAP
    w = ustrip(u"Å", solvent_layer_width)
    Lz = n_leaflets * (zsolvent + w)

    lipid_placements = _LipidPlacement[]
    solvent_slabs = Tuple{Float64,Float64}[]
    if type == :bilayer
        for upper in (true, false)
            s = upper ? 1.0 : -1.0
            for i in 1:nlip
                lipid_counts[i] == 0 && continue
                push!(lipid_placements, _LipidPlacement(
                    lipid_pdbfiles[i], lipid_counts[i], lipid_head[i], lipid_tail[i],
                    s * head_planes[i], s * _CHAIN_END_PLANE, upper, nothing,
                ))
            end
        end
        push!(solvent_slabs, (zsolvent, Lz / 2), (-Lz / 2, -zsolvent))
    else
        bottom = -Lz / 2
        for i in 1:nlip
            lipid_counts[i] == 0 && continue
            push!(lipid_placements, _LipidPlacement(
                lipid_pdbfiles[i], lipid_counts[i], lipid_head[i], lipid_tail[i],
                bottom + head_planes[i], bottom + _CHAIN_END_PLANE, true, bottom,
            ))
        end
        push!(solvent_slabs, (bottom + zsolvent, Lz / 2))
    end

    # Solvent: one slab per leaflet, all with the same composition and counts
    mean_M = sum(solvent_weight[k] * solvent_molar_mass[k] for k in 1:nsolv)
    slab_volume = Lxy^2 * solvent_layer_width
    total_mass = uconvert(u"g", solvent_density * slab_volume)
    total_moles = total_mass / mean_M
    solvent_counts = [round(Int, ustrip(solvent_weight[k] * total_moles * Unitful.Na)) for k in 1:nsolv]

    solvent_placements = _SolventPlacement[]
    for (zlo, zhi) in solvent_slabs
        for k in 1:nsolv
            solvent_counts[k] == 0 && continue
            push!(solvent_placements, _SolventPlacement(solvent_pdbfiles[k], solvent_counts[k], zlo, zhi))
        end
    end

    unitcell = Matrix{Float64}(Diagonal([ustrip(u"Å", Lxy), ustrip(u"Å", Lxy), Lz]))

    lipid_lines = join(
        ("    $(lipid_names[i]): $(lipid_counts[i]) per leaflet, " *
         "head-to-tail length = $(round(L[i]; digits=2)) Å" for i in 1:nlip),
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
        Area per lipid = $(round(u"Å^2", area_per_lipid; digits=2))
        Box lateral area = $(round(u"Å^2", total_area_used; digits=2))
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
    input, output = expanduser(input), expanduser(output)
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
            head, tail = lp.upper ? ("above", "below") : ("below", "above")
            floor_line = isnothing(lp.floor) ? "" : "    above plane 0. 0. 1. $(lp.floor)\n"
            println(io,
                """
                structure $(lp.pdbfile)
                    number $(lp.number)
                $(floor_line)    atoms $(lp.head_idx)
                        $head plane 0. 0. 1. $(lp.head_plane)
                    end atoms
                    atoms $(join(lp.tail_idx, " "))
                        $tail plane 0. 0. 1. $(lp.tail_plane)
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

Builds and packs a lipid membrane system, without writing an input file.

It is equivalent to calling
[`write_packmol_input`](@ref write_packmol_input(::Membrane)) followed by `packmol` on the
resulting file.

`output` behaves as in `write_packmol_input`. Any other keyword (`nloop`, `iprint`,
`seed`, ...) is forwarded to the packing engine — see `packmol(::PackmolSystem)`.

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
    z = [0.0, 0.0, 1.0]
    for lp in lipid_placements
        head_c = lp.upper ? AbovePlane(z, lp.head_plane) : BelowPlane(z, lp.head_plane)
        tail_c = lp.upper ? BelowPlane(z, lp.tail_plane) : AbovePlane(z, lp.tail_plane)
        # (a vector of structure types for a file of several conformations)
        st = structure_type(lp.pdbfile;
            number=lp.number,
            constraints=isnothing(lp.floor) ? AnyConstraint{Float64}[] : [AbovePlane(z, lp.floor)],
            atom_constraints=[[lp.head_idx] => [head_c], lp.tail_idx => [tail_c]],
        )
        st isa AbstractVector ? append!(structure_types, st) : push!(structure_types, st)
    end
    for sp in solvent_placements
        push!(structure_types, structure_type(sp.pdbfile; number=sp.number,
            constraints=[AbovePlane(z, sp.zlo), BelowPlane(z, sp.zhi)],
        ))
    end
    packmol_system = PackmolSystem(structure_types;
        output, tolerance=2.0, add_box_sides=true, seed=-1, _recipe_unitcell(unitcell)...,
    )
    return packmol(packmol_system; kwargs...)
end
