#
# The polymers in the Julia API: `polymer_structure`, the equivalent of a structure block
# with `repeat` blocks, and `monomer`, a monomer of its sequence.
#

"""
    monomer(natoms; links = [])

A monomer of the sequence of a polymer (see [`polymer_structure`](@ref)): its number of
atoms, in the order of the structure of the polymer, and its links to the next monomer,
each `(a, b)` (atom `a` of this monomer to atom `b` of the next one, at their distance in
the structure of the polymer), `(a, b, d0)`, `(a, b, d0, k)` or `(a, b, d0, k, w)`: the
distance, the constant (default $(_default_link_constant)) and the half-width of a flat
bottom (no penalty within `d0 ± w`; default 0) of the link, whose penalty is
`k (|d - d0| - w)²`.
"""
function monomer(natoms::Integer; links = Tuple[])
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
    return MonomerSpec{Float64}(Int(natoms), parsed)
end

# A polymer of the Julia API: the structure types of its conformations, and its sequence
struct PolymerStructure{D,T}
    rigid_types::Vector{StructureType{D,T}}
    spec::PolymerSpec{T}
end

"""
    polymer_structure(filename; number, sequence, cyclic = false, kwargs...)

The polymers of the structure file `filename` (a file with several models, or a vector of
files, for several conformations), as a structure block with `repeat` blocks: `number`
polymers, of the `sequence` of monomers given as `repeats => monomers` pairs, for
example `[20 => [monomer(18; links = [(17, 1), (18, 1), (18, 2)])]]` for 20 repeats of a
monomer of 18 atoms. With `cyclic = true`, the last monomer is linked to the first one.
The other keyword arguments are those of [`structure_type`](@ref) (`constraints`,
`atom_constraints`, `radius`, `radscale`, ...). Give it to `PackmolSystem` with the other
structure types, as `PackmolSystem([polymers, water]; ...)`.

Each monomer is packed as a rigid molecule, held to the next one by its links. The
polymers are first packed as rigid molecules, in their conformations, and then their
monomers are released.
"""
function polymer_structure(filename::Union{String,AbstractVector{String}};
    number::Int, sequence::AbstractVector{<:Pair}, cyclic::Bool = false, D::Int = 3, T::DataType = Float64, kwargs...,
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
            push!(monomers, MonomerSpec{T}(m.natoms, [PolymerLink{T}(l.a, l.b, l.d0, l.k, l.w) for l in m.links]))
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
    return PolymerStructure{D,T}(rigid_types, PolymerSpec{T}(monomers, cyclic))
end
