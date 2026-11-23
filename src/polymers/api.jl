#
# The flexible molecules in the Julia API: `flexible_structure_type`, the equivalent of a
# structure block with `fragment` lines (alone or in `repeat` blocks) or `rotatable`
# lines, and `fragment`, a fragment of its sequence.
#

"""
    fragment(natoms; links = [], rotatable = [])

A rigid fragment of a flexible molecule. It is given by its number of atoms (the next
`natoms` atoms of the structure of the molecule), and by its rotatable bonds or links to
the next fragment, in the sequence of fragments of [`flexible_structure_type`](@ref). The
fragments of a polymer are its monomers.

`rotatable` are bonds to the next fragment, each `(a, b)` (atom `a` of this fragment
bonded to atom `b` of the next one): each stands for the links that keep the length and
the angles of the bond (the bond, and the bonds of `a` and `b` to their neighbors in their
fragments), such that the fragments turn about it. They require the bonds of the
molecule: the CONECT records of its file, or `bond_distance` in `flexible_structure_type`.
The atom indices are local to each fragment, from 1 (the position of the atom among the
atoms of the fragment).

`links` are general links to the next fragment, each `(a, b)` (atom `a` of this fragment
to atom `b` of the next one, at their distance in the structure of the molecule),
`(a, b, d0)`, `(a, b, d0, k)` or `(a, b, d0, k, w)`: the distance, the constant (default
$(_default_link_constant)) and the half-width of a flat bottom (no penalty within
`d0 ± w`; default 0) of the link, whose penalty is `k (|d - d0| - w)²`. A single link is a
joint free to turn in any direction (the angles are not kept), and links need not join
bonded atoms.
"""
function fragment(natoms::Integer; links = Tuple[], rotatable = Tuple[])
    natoms >= 1 || throw(ArgumentError("""\n
        A fragment must have at least one atom.

    """))
    parsed = PolymerLink{Float64}[]
    for l in links
        (2 <= length(l) <= 5) || throw(ArgumentError("""\n
            A link is (a, b), (a, b, d0), (a, b, d0, k) or (a, b, d0, k, w): got $l.

        """))
        a, b = Int(l[1]), Int(l[2])
        (1 <= a <= natoms) || throw(ArgumentError("""\n
            link $l: atom $a is not an atom of the fragment, of $natoms atoms.

        """))
        value(i) = length(l) >= i ? Float64(l[i]) : NaN
        push!(parsed, PolymerLink{Float64}(a, b, value(3), value(4), value(5)))
    end
    bonds = Tuple{Int,Int}[]
    for r in rotatable
        length(r) == 2 || throw(ArgumentError("""\n
            A rotatable bond is (a, b): got $r.

        """))
        a, b = Int(r[1]), Int(r[2])
        (1 <= a <= natoms) || throw(ArgumentError("""\n
            rotatable $r: atom $a is not an atom of the fragment, of $natoms atoms.

        """))
        push!(bonds, (a, b))
    end
    return FragmentSpec{Float64}(Int(natoms), parsed, bonds)
end

# A flexible molecule of the Julia API: the structure types of its conformations (as
# rigid molecules), and its fragments or rotatable bonds
struct FlexibleStructureType{D,T}
    rigid_types::Vector{StructureType{D,T}}
    spec::PolymerSpec{T}
end

"""
    flexible_structure_type(filename; number, fragments = nothing, rotatable = nothing, cyclic = false, bond_distance = nothing, kwargs...)

The structure type of `number` flexible molecules. The molecules of the structure file
`filename` (a file with several models, or a vector of files, for several conformations)
are made of rigid fragments that turn about the bonds between them (or are held by general
links), as in a structure block with `fragment` or `rotatable` lines. Give either:

- `fragments`: the sequence of the fragments of the molecule, in the order of its atoms
  (see [`fragment`](@ref)), each a fragment or `n => fragments`, `n` repeats of a fragment
  or of a vector of them (a `repeat` block). For example, `[20 => fragment(18; rotatable =
  [(18, 1)])]` for a polymer of 20 monomers of 18 atoms, bonded by atom 18 of each to atom 1
  of the next, or `[fragment(5; rotatable = [(5, 1)]), 20 => [...], fragment(7)]` with end
  groups. With `cyclic = true`, the last fragment is linked to the first one.
- `rotatable`: the rotatable bonds of the molecule, a vector of pairs of bonded atoms, as
  `[(21, 22), (36, 39)]`, indexed by their positions in the file, from 1, or `:all`: the
  bonds not in rings with another heavy atom bonded to each end, except those between two
  carbon or nitrogen atoms with three bonds each (double, conjugated and amide bonds; for
  all-atom models). The molecule is cut at these bonds into its fragments.

The rotatable bonds require the bonds of the molecule: those of the CONECT records of the
file or, with `bond_distance = d`, the pairs of atoms closer than `d` Å (in the first
conformation). The other keyword arguments are those of [`structure_type`](@ref)
(`constraints`, `atom_constraints`, `radius`, `radscale`, ...; not `fixed`,
`restart_from` or `constrain_rotation`). Give it to `PackmolSystem` with the other
structure types, as `PackmolSystem([chains, water]; ...)`.

The molecules are first packed as rigid molecules, in their conformations, and then their
fragments are released.
"""
function flexible_structure_type(filename::Union{String,AbstractVector{String}};
    number::Int,
    fragments::Union{Nothing,AbstractVector} = nothing,
    rotatable::Union{Nothing,Symbol,AbstractVector} = nothing,
    cyclic::Bool = false,
    bond_distance::Union{Nothing,Real} = nothing,
    D::Int = 3, T::DataType = Float64, kwargs...,
)
    isnothing(fragments) == isnothing(rotatable) && throw(ArgumentError("""\n
        flexible_structure_type: give either `fragments` (the sequence of the fragments of the
        molecule) or `rotatable` (its rotatable bonds).

    """))
    cyclic && isnothing(fragments) && throw(ArgumentError("""\n
        flexible_structure_type: `cyclic` requires `fragments`.

    """))
    types = structure_type(filename; number, D, T, kwargs...)
    rigid_types = types isa StructureType ? [types] : types
    bond_d = T(something(bond_distance, NaN))
    if !isnothing(rotatable)
        return FlexibleStructureType{D,T}(rigid_types, PolymerSpec{T}(FragmentSpec{T}[], false, _rotatable_pairs(rotatable), bond_d))
    end
    specs = FragmentSpec{T}[]
    for entry in fragments
        if entry isa FragmentSpec
            push!(specs, _fragment_copy(T, entry))
        elseif entry isa Pair
            nrepeat, block = entry
            (nrepeat isa Integer && nrepeat >= 1) || throw(ArgumentError("""\n
                The number of repeats must be an integer of at least 1, got $nrepeat.

            """))
            block isa FragmentSpec && (block = [block])
            (block isa AbstractVector && !isempty(block) && all(m -> m isa FragmentSpec, block)) || throw(ArgumentError("""\n
                A repeat is `n => fragment(...)` or `n => [fragment(...), ...]`.

            """))
            for _ in 1:nrepeat, m in block
                push!(specs, _fragment_copy(T, m))
            end
        else
            throw(ArgumentError("""\n
                The fragments are `fragment(...)` or `n => fragments` (n repeats), got $(typeof(entry)).

            """))
        end
    end
    isempty(specs) && throw(ArgumentError("""\n
        flexible_structure_type: no fragments given.

    """))
    natoms = sum(m.natoms for m in specs)
    natoms == first(rigid_types).natoms || throw(ArgumentError("""\n
        The fragments of $(basename(first(rigid_types).filename)) have $natoms atoms, and its
        structure has $(first(rigid_types).natoms).

    """))
    return FlexibleStructureType{D,T}(rigid_types, PolymerSpec{T}(specs, cyclic, Tuple{Int,Int}[], bond_d))
end

_fragment_copy(::Type{T}, m::FragmentSpec) where {T} =
    FragmentSpec{T}(m.natoms, [PolymerLink{T}(l.a, l.b, l.d0, l.k, l.w) for l in m.links], copy(m.rotatable))

# The rotatable bonds given in the Julia API: pairs of atoms, or `:all` (the sentinel (0, 0))
function _rotatable_pairs(rotatable::Symbol)
    rotatable == :all || throw(ArgumentError("""\n
        flexible_structure_type: rotatable is a vector of pairs of atoms, or :all.

    """))
    return [(0, 0)]
end
function _rotatable_pairs(rotatable::AbstractVector)
    isempty(rotatable) && throw(ArgumentError("""\n
        flexible_structure_type: no rotatable bonds given.

    """))
    return Tuple{Int,Int}[(Int(b[1]), Int(b[2])) for b in rotatable]
end
