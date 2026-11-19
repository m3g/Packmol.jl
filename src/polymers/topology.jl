#
# The topology of the polymers of a system: the links between atoms of the system, the
# atoms of the excluded pairs, and the polymer of each monomer molecule.
#

# The index of the first molecule and of the first atom of each structure type, minus one
function _type_offsets(sys::PackmolSystem)
    ntypes = length(sys.structure_types)
    mol_offset = zeros(Int, ntypes)
    atom_offset = zeros(Int, ntypes)
    imol = 0
    iat = 0
    for (t, st) in enumerate(sys.structure_types)
        mol_offset[t] = imol
        atom_offset[t] = iat
        imol += st.number_of_molecules
        iat += st.number_of_molecules * st.natoms
    end
    return mol_offset, atom_offset
end

# The polymer (index in `sys.polymers`), the conformation and the monomer of each
# structure type (0, 0, 0 if not a monomer)
function _type_monomers(sys::PackmolSystem)
    ntypes = length(sys.structure_types)
    polymer_of_type = zeros(Int, ntypes)
    conformation_of_type = zeros(Int, ntypes)
    monomer_of_type = zeros(Int, ntypes)
    for (ip, p) in enumerate(sys.polymers)
        for c in 1:_nconformations(p), m in 1:_nmonomers(p)
            t = _monomer_type(p, c, m)
            polymer_of_type[t] = ip
            conformation_of_type[t] = c
            monomer_of_type[t] = m
        end
    end
    return polymer_of_type, conformation_of_type, monomer_of_type
end

# The molecule of monomer `m` of copy `j` of the polymers of conformation `c`
_monomer_molecule(p::Polymer, mol_offset::Vector{Int}, c::Int, m::Int, j::Int) = mol_offset[_monomer_type(p, c, m)] + j

#
# Sets the topology of the polymers of the system (`sys.polymer_topology`)
#
function _set_polymer_topology!(sys::PackmolSystem{D,T}) where {D,T}
    if isempty(sys.polymers)
        sys.polymer_topology = PolymerTopology{T}()
        return sys
    end
    mol_offset, atom_offset = _type_offsets(sys)
    links = SystemLink{T}[]
    template_atom = zeros(Int32, length(sys.atoms))
    mol_instance = zeros(Int32, sys.nmols)
    instance_polymer = Int32[]
    instance_conformation = Int32[]
    for (ip, p) in enumerate(sys.polymers)
        nmon = _nmonomers(p)
        for (c, j) in _polymer_copies(p)
            push!(instance_polymer, ip)
            push!(instance_conformation, c)
            instance = Int32(length(instance_polymer))
            for m in 1:nmon
                t = _monomer_type(p, c, m)
                natm = sys.structure_types[t].natoms
                imol = mol_offset[t] + j
                mol_instance[imol] = instance
                iat0 = atom_offset[t] + (j - 1) * natm
                for (k, i) in enumerate(p.monomer_atoms[m])
                    p.excluded_atom[c][i] && (template_atom[iat0+k] = Int32(i))
                end
            end
            # the links between the fragments of this copy
            for l in p.links[c]
                ta = _monomer_type(p, c, l.fa)
                tb = _monomer_type(p, c, l.fb)
                ia = atom_offset[ta] + (j - 1) * sys.structure_types[ta].natoms + l.a
                ib = atom_offset[tb] + (j - 1) * sys.structure_types[tb].natoms + l.b
                push!(links, SystemLink{T}(ia, ib, mol_offset[ta] + j, mol_offset[tb] + j, l.d0, l.k, l.w))
            end
        end
    end
    sys.polymer_topology = PolymerTopology{T}(links, template_atom, mol_instance, instance_polymer, instance_conformation)
    return sys
end

# Whether molecule `imol` is a monomer of a polymer
_is_monomer(topo::PolymerTopology, imol::Integer) = !isempty(topo.mol_instance) && topo.mol_instance[imol] != 0

#
# Whether the pair of atoms `i`, `j` (of different molecules) is excluded from the
# distances: atoms of the same polymer, close in its structure
#
@inline function _excluded_pair(sys::PackmolSystem, i::Int, j::Int, imol::Int, jmol::Int)
    topo = sys.polymer_topology
    isempty(topo.template_atom) && return false
    ti = topo.template_atom[i]
    ti == 0 && return false
    tj = topo.template_atom[j]
    tj == 0 && return false
    q = topo.mol_instance[imol]
    q == topo.mol_instance[jmol] || return false
    exclusions = sys.polymers[topo.instance_polymer[q]].exclusions[topo.instance_conformation[q]]
    return (min(ti, tj), max(ti, tj)) in exclusions
end

# The periodic image of the vector `v` closest to the origin (`v` reduced to the cell
# first: `_nearest_periodic_image` searches the neighboring images only, and the vectors
# between the monomers of an unwrapped chain can be several cells long)
_min_image(v::SVector{3,T}, unitcell::AbstractMatrix) where {T} =
    _nearest_periodic_image(wrap_to_center(v, unitcell, zero(v)), unitcell, zero(v))
_min_image(v::SVector{D,T}, unitcell::AbstractMatrix) where {D,T} = wrap_to_center(v, unitcell, zero(v))
