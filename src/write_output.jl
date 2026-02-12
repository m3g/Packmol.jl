"""
    write_output(packmol_system::PackmolSystem)

Write packed coordinates to the output file specified in `packmol_system.output_file`.
"""
function write_output(packmol_system::PackmolSystem{D,T}; output_file=packmol_system.output_file) where {D,T}
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    output_atoms = Atom[]
    iat = 0
    imol = 0
    atom_serial = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            for j in 1:st.natoms
                iat += 1
                atom_serial += 1
                atom = copy(st.atoms[j])
                pos = atom_positions[iat]
                atom.index = atom_serial
                atom.x = pos[1]
                atom.y = pos[2]
                if D == 3
                    atom.z = pos[3]
                end
                atom.resnum = imol
                push!(output_atoms, atom)
            end
        end
    end

    write_pdb(output_file, output_atoms)
    return nothing
end


