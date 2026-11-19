#
# The types of polymers (structure blocks with `repeat` blocks): the monomers given in
# the input, the links between consecutive monomers, and the topology of the links and
# of the excluded pairs of atoms in a system.
#
# A polymer is packed as a chain of rigid monomers: each monomer is a structure type
# (one molecule per polymer), consecutive monomers held together by harmonic links
# between their atoms, `k (d - d0)²`. The pairs of atoms of different monomers closer
# than the sum of their radii in the structure of the polymer (the bonded atoms of a
# junction, and their neighbors) are excluded from the distances.
#

# A link between atom `a` of a monomer and atom `b` of the next one, with its distance,
# constant and the half-width of its flat bottom (`NaN`: not given, the distance of the
# structure, the default constant, and no flat bottom)
struct PolymerLink{T}
    a::Int
    b::Int
    d0::T
    k::T
    w::T
end

# A monomer of the sequence of a polymer: its number of atoms, its links to the next, and
# its rotatable bonds to the next (atom of this monomer, atom of the next), each standing
# for the links of the bond (built from the bonds of the polymer)
struct MonomerSpec{T}
    natoms::Int
    links::Vector{PolymerLink{T}}
    rotatable::Vector{Tuple{Int,Int}}
end
MonomerSpec{T}(natoms::Int, links::Vector{PolymerLink{T}}) where {T} = MonomerSpec{T}(natoms, links, Tuple{Int,Int}[])

# The flexibility given in a structure block: the sequence of monomers of a polymer, or
# the rotatable bonds of a molecule, with the largest distance of bonded atoms if the
# bonds are computed from the distances (`NaN`: from the CONECT records); both empty:
# a rigid molecule
struct PolymerSpec{T}
    monomers::Vector{MonomerSpec{T}}
    cyclic::Bool
    rotatable::Vector{Tuple{Int,Int}}
    bond_distance::T
end
PolymerSpec{T}() where {T} = PolymerSpec{T}(MonomerSpec{T}[], false, Tuple{Int,Int}[], T(NaN))

# A link between atom `a` of fragment `fa` and atom `b` of fragment `fb` of a flexible
# molecule (indices of the atoms in their fragments), with its distance, constant and
# half-width of its flat bottom
struct FragmentLink{T}
    fa::Int
    a::Int
    fb::Int
    b::Int
    d0::T
    k::T
    w::T
end

#
# A flexible molecule of a system (a polymer, or a molecule with rotatable bonds): the
# structure types of its conformations as rigid molecules (`rigid_types`, for the initial
# placement), the index of the structure type of its first fragment, the atoms of each
# fragment in the structure of the molecule (the monomers of a polymer, in order), the
# links between the fragments in each conformation (`links[c]`), whether the polymer is
# cyclic (for its summary), and the pairs of atoms of the structure excluded from the
# distances in each conformation (`exclusions[c]`, `i < j`), with the atoms in any of
# them (`excluded_atom[c]`).
#
# The fragments are consecutive structure types, conformation after conformation: the
# fragments of a molecule have the geometry of its conformation (fragments of a
# conformation aligned to the atoms of another would not fit, and the links would be
# strained). The names keep the words of the polymers: a monomer is a fragment.
#
struct Polymer{D,T}
    rigid_types::Vector{StructureType{D,T}}
    first_type::Int
    monomer_atoms::Vector{Vector{Int}}
    links::Vector{Vector{FragmentLink{T}}}
    cyclic::Bool
    exclusions::Vector{Set{Tuple{Int32,Int32}}}
    excluded_atom::Vector{Vector{Bool}}
end

_nmonomers(p::Polymer) = length(p.monomer_atoms)
_nconformations(p::Polymer) = length(p.rigid_types)
_npolymers(p::Polymer) = sum(st.number_of_molecules for st in p.rigid_types; init=0)
_ntypes(p::Polymer) = _nconformations(p) * _nmonomers(p)
# The structure type of monomer `m` of the polymers of conformation `c`
_monomer_type(p::Polymer, c::Int, m::Int) = p.first_type + (c - 1) * _nmonomers(p) + m - 1
# The copies of a polymer, (conformation, copy of the conformation), in the order of the
# rigid molecules
function _polymer_copies(p::Polymer)
    copies = Tuple{Int,Int}[]
    for (c, st) in enumerate(p.rigid_types), j in 1:st.number_of_molecules
        push!(copies, (c, j))
    end
    return copies
end

# A link between two atoms of the system, of the molecules `ma` and `mb`: no penalty
# for distances within `d0 ± w` (a flat bottom), `k (|d - d0| - w)²` outside
struct SystemLink{T}
    ia::Int
    ib::Int
    ma::Int
    mb::Int
    d0::T
    k::T
    w::T
end

#
# The links and the excluded pairs of the polymers of a system (empty without polymers):
# `template_atom[i]` is the index, in the structure of its polymer, of atom `i` of the
# system if it is in an excluded pair (0 otherwise); `mol_instance[imol]` is the polymer
# (counted over all the polymers of the system) of molecule `imol` (0 if it is not a
# monomer), `instance_polymer[q]` the index of polymer `q` in `polymers`, and
# `instance_conformation[q]` its conformation.
#
struct PolymerTopology{T}
    links::Vector{SystemLink{T}}
    template_atom::Vector{Int32}
    mol_instance::Vector{Int32}
    instance_polymer::Vector{Int32}
    instance_conformation::Vector{Int32}
end
PolymerTopology{T}() where {T} = PolymerTopology{T}(SystemLink{T}[], Int32[], Int32[], Int32[], Int32[])

# One-line summary of a polymer (in the `show` of a `PackmolSystem`)
function _summary(p::Polymer)
    st = first(p.rigid_types)
    nconf = _nconformations(p)
    return string(
        "\"", basename(st.filename), "\" - flexible: ", _nmonomers(p), " fragments",
        p.cyclic ? " (cyclic)" : "", nconf > 1 ? " - $nconf conformations" : "",
        " - n_atoms = ", st.natoms, " - n_mols = ", _npolymers(p), " - n_constraints = ", length(st.constraints),
    )
end
