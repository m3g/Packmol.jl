#
# The system to be packed (`PackmolSystem`): its structure types, atoms, variables and
# options; its constructors, `show`, and the scale of the radii of each structure type.
#

@kwdef mutable struct PackmolSystem{D,T}
    filetype::String = "pdb"
    input_file::String
    output_file::String
    tolerance::T = 2.0
    # Loose-start radius scale factor for the packing objective (Fortran Packmol's
    # `discale`): atom radii are inflated by this factor at the start of packing
    # (loop 0), then decayed exponentially toward 1.0 by loop `radscale_decay_loops`, the
    # same for every atom regardless of movebad! relocations (see
    # `_radscale_at` and the outer packing loop in packing/packing_loop.jl). `nothing` (the
    # default) means a scale for each structure type (`_radscale`): 5.0 for the free
    # molecules of more than 20 atoms, which otherwise interpenetrate (lipids, liquid
    # crystals: 2 to 6 times faster, and 18,000 lipids of an exosome that did not pack
    # at 1.0 packed in 8 loops), unless restarted (`restart_from`), 1.0 for the smaller ones (water, ions, small solvents
    # pack in a few loops at 1.0, and 4 to 8 times slower inflated), and the largest scale
    # of the free molecules for the fixed ones (lipids inflated next to a fixed protein at
    # its radii were pushed into it). A `radscale` keyword in a structure block sets the
    # scale of that type.
    radscale::Union{Nothing,T} = nothing
    # Loops for the radii to decay from `radscale` to 1.0 (input keyword
    # `radscale_decay_loops`): `nothing` means 20.
    radscale_decay_loops::Union{Nothing,Int} = nothing
    # Minimum time, in seconds, between two writes of the output file during the
    # packing (input keyword `writeout`, as in Fortran Packmol): `nothing` means 60 s.
    writeout::Union{Nothing,T} = nothing
    structure_types::Vector{StructureType{D,T}} = StructureType{D,T}[]
    tolerance_precision::T = 1e-1
    constraint_precision::T = 1e-2
    # Packing loops and optimizer iterations per loop (input keywords `nloop`
    # and `maxit`): `nothing` means unset, and the keyword arguments of
    # `packmol(::PackmolSystem)` take precedence over these.
    nloop::Union{Nothing,Int} = nothing
    max_iter::Union{Nothing,Int} = nothing
    # Maximum fraction of the molecules of each type moved by `movebad!` per
    # loop (input keyword `movefrac`): `nothing` means unset, as above.
    movefrac::Union{Nothing,T} = nothing
    # Maximum number of threads used (input keyword `nthreads`): `nothing` means all.
    nthreads::Union{Nothing,Int} = nothing
    max_random_init::Int=20
    adjust_constraints_on_init::Bool=true
    add_amber_ter::Bool = false
    amber_ter_preserve::Bool = false
    add_box_sides::Bool = false
    # Write the connectivity of the structure files (CONECT records) to the output;
    # disabled with `ignore_conect` (or `connect no`).
    connect::Bool = true
    seed::Int = 1234567
    avoid_overlap::Bool = true
    # Print the report of the packing (and its warnings) to the screen or log; `false`
    # (input keyword `screen_log no`) silences the run
    screen_log::Bool = true
    # Write the initial approximation to the output file and stop, without packing (input
    # keyword `check_initial_placement`, or Fortran's `check`)
    check_initial_placement::Bool = false
    # Pack even if the initial approximation does not satisfy the constraints (input
    # keyword `continue_failed_initial_packing`); otherwise `packmol` stops there.
    continue_failed_initial_packing::Bool = false
    # Whole-system restart (Fortran's `restart_from` with no structure type,
    # i.e. index 0): `nothing` means unset. When
    # `restart_from` is set, packmol() skips the whole initial-approximation
    # pipeline entirely and reads every molecule's position from this
    # source instead — the point of a whole-system restart is to skip that
    # (potentially expensive, for a large system) phase, not just to get
    # the same final starting point some other way.
    # `restart_from` accepts either a coordinate file (PDB, XYZ or mmCIF, by
    # extension), whose atoms are rigid-body-aligned against each structure
    # type's own template to recover (cm, angles) (see `_align_molecule` in
    # write_output.jl), or, Julia-API only, a `Vector{<:Atom}` already in
    # memory. There is no `restart_to`: the output file itself can be given to
    # `restart_from`. See also `StructureType.restart_from` for the
    # per-structure-type equivalent.
    restart_from::Union{Nothing,String,Vector{Atom{Nothing}}} = nothing
    # Unit cell for periodic boundary conditions (nothing = no PBC)
    unitcell::Union{Nothing, Matrix{T}} = nothing
    # Reference center for PBC wrapping (constraints evaluated relative to this point)
    unitcell_center::Union{Nothing, SVector{D,T}} = nothing
    # Which fundamental-domain shape `get_atoms`/`write_output` wrap molecule
    # positions into, when `unitcell` is set: `:triclinic` (the ordinary
    # skewed-parallelepiped shape, e.g. what CRYST1/PDB output implies),
    # `:dodecahedral` (the rhombic-dodecahedron/Wigner-Seitz "compact"
    # shape), or `:octahedral` (the truncated-octahedron/Wigner-Seitz
    # "compact" shape). Set via `triclinic_to_dodecahedral(::PackmolSystem)`/
    # `dodecahedral_to_triclinic(::PackmolSystem)` (or their `octahedral`
    # counterparts) rather than directly.
    periodic_boundary_style::Symbol = :triclinic
    # Whether the cell is dodecahedral or octahedral (see `_is_wigner_seitz_cell`): the
    # constraints are then evaluated, and the molecules confined, in the rhombic
    # dodecahedron or truncated octahedron centered at `unitcell_center`, instead of in
    # the parallelepiped of the cell vectors. Set with the confinement constraints.
    wigner_seitz_cell::Bool = false
    # Internal data for the optimization
    nmols::Int = 0
    atoms::Vector{AtomData{T}} = AtomData{T}[]
    molecule_positions::Vector{MoleculePosition{D,T}} = MoleculePosition{D,T}[]
    # Outcome of the last `packmol(::PackmolSystem)` call on this system:
    # `:not_packed` (never packed, e.g. freshly built, or only `check`-mode
    # placement), `:packing_ready` (converged within `tolerance_precision`/
    # `constraint_precision`), or `:failed_packing` (ran out of `nloop`
    # loops short of that). `packmol` returns this same `PackmolSystem`
    # (mutated in place) rather than a bare `Bool`, precisely so this field
    # is how callers check the outcome — see `packmol(::PackmolSystem)`.
    status::Symbol = :not_packed
end

#
# Build the per-atom (`AtomData`) and per-molecule (`MoleculePosition`) arrays
# from a vector of `StructureType`s: one is a flat list of every atom across
# every molecule of every structure type, tagging it with its molecule index,
# structure type index, radius, and constraint indices; the other is a
# placeholder (zeroed) position per molecule, later filled in by
# `initialize_molecules!`. Shared by `read_packmol_input` and the
# `PackmolSystem` outer constructor below.
#
function _atoms_and_molecule_positions(structure_types::Vector{StructureType{D,T}}) where {D,T}
    mol_index = 0
    atoms = AtomData{T}[]
    molecule_positions = MoleculePosition{D,T}[]
    for (itype, st) in enumerate(structure_types)
        for _ in 1:st.number_of_molecules
            mol_index += 1
            push!(molecule_positions, MoleculePosition(zeros(SVector{D,T}), zeros(SVector{D,T})))
            for iatom in 1:st.natoms
                push!(atoms, AtomData(mol_index, itype, st.radii[iatom], st.atom_constraints[iatom]))
            end
        end
    end
    return atoms, molecule_positions, mol_index
end

#
# Append one constraint to a structure type's own constraint list, wiring it
# into every atom's `atom_constraints` (not just the atoms present when the
# structure type was originally built) — including atoms whose own
# `atoms ... end atoms` block already narrowed them to a custom subset of
# constraints, since this one is a whole-system requirement, not a per-atom
# feature toggle.
#
function _push_constraint!(st::StructureType{D,T}, c::AnyConstraint{T}) where {D,T}
    push!(st.constraints, c)
    idx = length(st.constraints)
    for ac in st.atom_constraints
        push!(ac, idx)
    end
    return st
end

"""
    PackmolSystem(structure_types::Vector{<:StructureType}; output::String, kwargs...)

Build a `PackmolSystem` directly from a vector of `StructureType`s (see
[`structure_type`](@ref)), without going through an input file. `output` sets
the output file path; any other `PackmolSystem` field (`tolerance`, `seed`,
`unitcell`, ...) can be given as a keyword argument.
"""
function PackmolSystem(
    structure_types::Vector{StructureType{D,T}};
    output::String,
    input_file::String = "",
    kwargs...,
) where {D,T}
    unitcell_raw = get(kwargs, :unitcell, nothing)
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(structure_types)
    kw = values(kwargs)
    if haskey(kw, :restart_from)
        kw = merge(kw, (; restart_from=_restart_source(kw.restart_from)))
    end
    # A unit cell without a center is centered at the origin (as `unitcell` in input files)
    if !isnothing(unitcell_raw) && !haskey(kw, :unitcell_center)
        kw = merge(kw, (; unitcell_center=zero(SVector{D,T})))
    end
    if !isnothing(unitcell_raw)
        kw = merge(kw, (; wigner_seitz_cell=_is_wigner_seitz_cell(Matrix{T}(unitcell_raw))))
    end
    return PackmolSystem{D,T}(;
        structure_types, atoms, molecule_positions, nmols,
        input_file, output_file=output, kw...,
    )
end

# A vector of structure types and of vectors of them (the conformations of a molecule,
# from `structure_type`), as `[lipids, water]`
function PackmolSystem(structure_types::AbstractVector; kwargs...)
    sts = StructureType[]
    for s in structure_types
        if s isa StructureType
            push!(sts, s)
        elseif s isa AbstractVector{<:StructureType}
            append!(sts, s)
        else
            throw(ArgumentError("""\n
                PackmolSystem: expected structure types (see `structure_type`), got a $(typeof(s)).

            """))
        end
    end
    isempty(sts) && throw(ArgumentError("PackmolSystem: no structure types given."))
    sts = [s for s in sts] # narrowed to Vector{StructureType{D,T}}
    isconcretetype(eltype(sts)) || throw(ArgumentError("PackmolSystem: structure types of different dimensions or precisions."))
    return PackmolSystem(sts; kwargs...)
end

function _indent(s::AbstractString; n=4)
    indented_str = IOBuffer()
    for line in eachline(IOBuffer(s))
        println(indented_str, repeat(" ", n) * line)
    end
    return String(take!(indented_str))
end

function Base.show(io::IO, ::MIME"text/plain", sys::PackmolSystem{D,T}) where {D,T}
    printstyled(io, "PackmolSystem{$D,$T}"; bold=true, color=:blue)
    if length(sys.input_file) > 0
        printstyled(io, " - read from: $(basename(sys.input_file))"; color=:blue)
    end
    println(io)
    groups = _molecule_groups(sys.structure_types)
    printstyled(io, _indent("Number of types of structures: $(length(groups))\n"); bold=true)
    printstyled(io, _indent("Total number of molecules: $(sys.nmols)\n"); bold=true)
    printstyled(io, _indent("Structure types:\n"); bold=true)
    for (i, group) in enumerate(groups)
        print(io, "    $i. ", _summary(@view(sys.structure_types[group])))
        i < length(groups) && println(io)
    end
end

#
# Scale of the radii of the atoms of a structure type at the start of the packing: the
# `radscale` of the structure block, or the global one, or, if neither is set, 5.0 for
# the free molecules of more than 20 atoms, 1.0 for the other free molecules and for the
# restarted ones (`restart_from`, of the system or of the type: inflated radii undo the
# arrangement of a restart, as a packed configuration or one built by another tool), and,
# for the fixed molecules, the largest scale of the free ones (see the `radscale` field of
# `PackmolSystem`). The fixed molecules are inflated as the free ones: with lipids
# inflated 5 times next to a fixed protein at its radii, the targets of the lipid-protein
# pairs (5 + 1 radii) were much shorter than those of the lipid-lipid pairs (5 + 5), and
# the lipids, crowded at their inflated radii, were pushed deep into the protein.
#
const _radscale_default = 5.0
const _radscale_min_natoms = 21
function _radscale(sys::PackmolSystem{D,T}, st::StructureType) where {D,T}
    default = if st.fixed.fixed
        maximum((_radscale_free(sys, s) for s in sys.structure_types if !s.fixed.fixed); init=one(T))
    else
        _radscale_free(sys, st)
    end
    return T(something(st.radscale, sys.radscale, default))
end
function _radscale_free(sys::PackmolSystem{D,T}, st::StructureType) where {D,T}
    restarted = !isnothing(sys.restart_from) || !isnothing(st.restart_from)
    default = (st.natoms >= _radscale_min_natoms && !restarted) ? T(_radscale_default) : one(T)
    return T(something(st.radscale, sys.radscale, default))
end
# Scale of the distance to the fixed atoms in the initial placement and in movebad!:
# the global radscale, if set (1.0 otherwise: free molecules are placed and moved at
# the tolerance distance from the fixed ones)
_global_radscale(sys::PackmolSystem{D,T}) where {D,T} = T(something(sys.radscale, one(T)))

#
