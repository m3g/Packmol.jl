#
# The polymers in the Julia API: `polymer_structure`, the equivalent of a structure block
# with `repeat` blocks, and `monomer`, a monomer of its sequence.
#

"""
    monomer(natoms; links = [], rotatable = [])

A monomer of the sequence of a polymer (see [`polymer_structure`](@ref)): its number of
atoms, in the order of the structure of the polymer, and its links to the next monomer,
each `(a, b)` (atom `a` of this monomer to atom `b` of the next one, at their distance in
the structure of the polymer), `(a, b, d0)`, `(a, b, d0, k)` or `(a, b, d0, k, w)`: the
distance, the constant (default $(_default_link_constant)) and the half-width of a flat
bottom (no penalty within `d0 ± w`; default 0) of the link, whose penalty is
`k (|d - d0| - w)²`. The atom indices are local to each monomer, from 1 (the position of
the atom among the atoms of the monomer).

`rotatable` are bonds to the next monomer, each `(a, b)` (atom `a` of this monomer bonded
to atom `b` of the next one), in place of links: each stands for the links that keep the
length and the angles of the bond (the bond, and the bonds of `a` and `b` to their
neighbors in their monomers), such that the monomers turn about it. They require the
bonds of the polymer: the CONECT records of its file, or `bond_distance` in
`polymer_structure`.
"""
function monomer(natoms::Integer; links = Tuple[], rotatable = Tuple[])
    natoms >= 1 || throw(ArgumentError("""\n
        A monomer must have at least one atom.

    """))
    parsed = PolymerLink{Float64}[]
    for l in links
        (2 <= length(l) <= 5) || throw(ArgumentError("""\n
            A link is (a, b), (a, b, d0), (a, b, d0, k) or (a, b, d0, k, w): got $l.

        """))
        a, b = Int(l[1]), Int(l[2])
        (1 <= a <= natoms) || throw(ArgumentError("""\n
            link $l: atom $a is not an atom of the monomer, of $natoms atoms.

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
            rotatable $r: atom $a is not an atom of the monomer, of $natoms atoms.

        """))
        push!(bonds, (a, b))
    end
    return MonomerSpec{Float64}(Int(natoms), parsed, bonds)
end

# A polymer of the Julia API: the structure types of its conformations, and its sequence
struct PolymerStructure{D,T}
    rigid_types::Vector{StructureType{D,T}}
    spec::PolymerSpec{T}
end

"""
    polymer_structure(filename; number, sequence, cyclic = false, bond_distance = nothing, kwargs...)

The polymers of the structure file `filename` (a file with several models, or a vector of
files, for several conformations), as a structure block with `repeat` blocks: `number`
polymers, of the `sequence` of monomers given as `repeats => monomers` pairs, for
example `[20 => [monomer(18; links = [(17, 1), (18, 1), (18, 2)])]]` for 20 repeats of a
monomer of 18 atoms. With `cyclic = true`, the last monomer is linked to the first one.
The rotatable bonds of the monomers (`monomer(natoms; rotatable)`) require the bonds of the
polymer: those of the CONECT records of the file or, with `bond_distance = d`, the pairs of
atoms closer than `d` Å (in the first conformation).
The other keyword arguments are those of [`structure_type`](@ref) (`constraints`,
`atom_constraints`, `radius`, `radscale`, ...). Give it to `PackmolSystem` with the other
structure types, as `PackmolSystem([polymers, water]; ...)`.

Each monomer is packed as a rigid molecule, held to the next one by its links. The
polymers are first packed as rigid molecules, in their conformations, and then their
monomers are released.
"""
function polymer_structure(filename::Union{String,AbstractVector{String}};
    number::Int, sequence::AbstractVector{<:Pair}, cyclic::Bool = false, bond_distance::Union{Nothing,Real} = nothing,
    D::Int = 3, T::DataType = Float64, kwargs...,
)
    types = structure_type(filename; number, D, T, kwargs...)
    rigid_types = types isa StructureType ? [types] : types
    monomers = MonomerSpec{T}[]
    for (nrepeat, block) in sequence
        nrepeat >= 1 || throw(ArgumentError("""\n
            The number of repeats must be at least 1, got $nrepeat.

        """))
        block isa MonomerSpec && (block = [block])
        for _ in 1:nrepeat, m in block
            push!(monomers, MonomerSpec{T}(m.natoms, [PolymerLink{T}(l.a, l.b, l.d0, l.k, l.w) for l in m.links], copy(m.rotatable)))
        end
    end
    isempty(monomers) && throw(ArgumentError("""\n
        polymer_structure: the sequence has no monomers.

    """))
    natoms = sum(m.natoms for m in monomers)
    natoms == first(rigid_types).natoms || throw(ArgumentError("""\n
        The monomers of the polymer $(basename(first(rigid_types).filename)) have $natoms atoms, and its
        structure has $(first(rigid_types).natoms).

    """))
    return PolymerStructure{D,T}(rigid_types, PolymerSpec{T}(monomers, cyclic, Tuple{Int,Int}[], T(something(bond_distance, NaN))))
end

"""
    flexible_structure(filename; number, rotatable, bond_distance = nothing, kwargs...)

The molecules of the structure file `filename` (a file with several models, or a vector
of files, for several conformations), flexible about the `rotatable` bonds, a vector of
pairs of bonded atoms, as `[(21, 22), (36, 39)]` (the `rotatable` lines of a structure
block), indexed by their positions in the file, from 1, or `:all` (`rotatable all`): the
bonds not in rings with another heavy atom bonded to each end, except those between two
carbon or nitrogen atoms with three bonds each (double, conjugated and amide bonds; for
all-atom models). The bonds of the molecule are those of the CONECT records of the file
or, with `bond_distance = d`, the pairs of atoms closer than `d` Å (in the first
conformation).
The other keyword arguments are those of [`structure_type`](@ref). Give it to
`PackmolSystem` with the other structure types, as `PackmolSystem([lipids, water]; ...)`.

The molecule is cut at the rotatable bonds into rigid fragments, kept together by links
that keep the length and the angles of each rotatable bond: the fragments turn about the
bonds only. The molecules are first packed as rigid molecules, in their conformations,
and then their fragments are released.
"""
function flexible_structure(filename::Union{String,AbstractVector{String}};
    number::Int, rotatable::Union{Symbol,AbstractVector}, bond_distance::Union{Nothing,Real} = nothing,
    D::Int = 3, T::DataType = Float64, kwargs...,
)
    types = structure_type(filename; number, D, T, kwargs...)
    rigid_types = types isa StructureType ? [types] : types
    if rotatable isa Symbol
        rotatable == :all || throw(ArgumentError("""\n
            flexible_structure: rotatable is a vector of pairs of atoms, or :all.

        """))
        pairs = [(0, 0)]
    else
        isempty(rotatable) && throw(ArgumentError("""\n
            flexible_structure: no rotatable bonds given.

        """))
        pairs = Tuple{Int,Int}[(Int(b[1]), Int(b[2])) for b in rotatable]
    end
    return PolymerStructure{D,T}(rigid_types, PolymerSpec{T}(MonomerSpec{T}[], false, pairs, T(something(bond_distance, NaN))))
end
