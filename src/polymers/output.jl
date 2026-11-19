#
# The output of the flexible molecules (polymers, molecules with rotatable bonds): each
# written whole, its atoms in the order of its structure, with a chain identifier per
# molecule (or the `chain` of its structure block),
# and the bonds of its structure; with PBC, the monomers of a polymer are put next to
# each other (each in the periodic image closest to the previous one), and the polymer
# wrapped into the cell by its center (the mean of the centers of its monomers), as rigid
# molecules are by their centers of mass.
#

#
# With PBC: shifts the atoms of the monomers of each polymer (`atom_positions`, unwrapped)
# such that the polymer is whole, its center wrapped into the cell by `wrap`
#
function _wrap_polymers!(
    atom_positions::Vector{SVector{D,T}}, sys::PackmolSystem{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}, wrap::F,
) where {D,T,F}
    mol_offset, atom_offset = _type_offsets(sys)
    for p in sys.polymers
        nmon = _nmonomers(p)
        targets = Vector{SVector{D,T}}(undef, nmon)
        for (c, j) in _polymer_copies(p)
            # the centers of the monomers, the chain made whole from the first one
            previous = zero(SVector{D,T})
            for m in 1:nmon
                cm = sys.molecule_positions[_monomer_molecule(p, mol_offset, c, m, j)].cm
                targets[m] = m == 1 ? cm : previous + _min_image(cm - previous, unitcell)
                previous = targets[m]
            end
            # the whole polymer wrapped by its center
            polymer_center = sum(targets) / nmon
            polymer_shift = wrap(polymer_center, unitcell, center) - polymer_center
            for m in 1:nmon
                t = _monomer_type(p, c, m)
                st = sys.structure_types[t]
                shift = targets[m] + polymer_shift - sys.molecule_positions[mol_offset[t]+j].cm
                iat0 = atom_offset[t] + (j - 1) * st.natoms
                for k in 1:st.natoms
                    atom_positions[iat0+k] += shift
                end
            end
        end
    end
    return atom_positions
end

#
# The atoms of the output (`atoms`, in the order of the structure types) with the
# polymers whole: the atoms of each copy of a polymer, monomer after monomer, at the
# place of its first monomer type, renumbered, with a chain per polymer
#
function _polymer_output_order(sys::PackmolSystem, atoms::Vector{Atom{Nothing}})
    isempty(sys.polymers) && return atoms
    polymer_of_type, _, _ = _type_monomers(sys)
    _, atom_offset = _type_offsets(sys)
    output = Atom{Nothing}[]
    ichain = 0
    for (t, st) in enumerate(sys.structure_types)
        ip = polymer_of_type[t]
        if ip == 0
            n = _ncopies(st) * st.natoms
            for i in atom_offset[t]+1:atom_offset[t]+n
                push!(output, atoms[i])
            end
        elseif t == sys.polymers[ip].first_type
            p = sys.polymers[ip]
            chain = first(p.rigid_types).chain
            # the atoms of each copy, from its fragments, in the order of the structure
            slots = Vector{Atom{Nothing}}(undef, first(p.rigid_types).natoms)
            for (c, j) in _polymer_copies(p)
                ichain += 1
                for m in 1:_nmonomers(p)
                    tm = _monomer_type(p, c, m)
                    natm = sys.structure_types[tm].natoms
                    i0 = atom_offset[tm] + (j - 1) * natm
                    for (k, i) in enumerate(p.monomer_atoms[m])
                        slots[i] = atoms[i0+k]
                    end
                end
                for atom in slots
                    atom.chain = string(isnothing(chain) ? _chainc(ichain) : chain)
                    push!(output, atom)
                end
            end
        end
    end
    for (i, atom) in enumerate(output)
        atom.index = i
    end
    return output
end

#
# The bonds of a polymer in the output: those of its structure, for each copy
#
function _polymer_connectivity!(connectivity::Vector{Vector{Int32}}, p::Polymer, offset::Int32)
    template = first(p.rigid_types).connectivity
    natoms = first(p.rigid_types).natoms
    for _ in 1:_npolymers(p)
        for iat in 1:natoms
            bonded = isnothing(template) || iat > length(template) ? Int32[] : template[iat] .+ offset
            push!(connectivity, bonded)
        end
        offset += Int32(natoms)
    end
    return offset
end
