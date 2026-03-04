#
# Pre-allocated scratch buffers for the initial-approximation loop.
# Created once (in `packmol`) and reused across all hot-path calls to avoid
# per-iteration heap allocation.
#
struct MemoryBuffers{D,T,V}
    atom_positions::Vector{SVector{D,T}}  # all-atom positions       (size: natoms)
    badness::Vector{T}                     # per-molecule badness     (size: nmols)
    free_atoms::Vector{SVector{D,T}}       # free-atom positions      (size: n_free_atoms)
    free_atom_mol::Vector{Int}             # free-atom → mol index    (size: n_free_atoms, stable)
    x::Vector{T}                           # optimizer variable       (size: nfree × 2D)
    vaux::V                                # SPGBox auxiliary storage (reused across loop iterations)
end

function MemoryBuffers(packmol_system::PackmolSystem{D,T}) where {D,T}
    natoms = length(packmol_system.atoms)
    nmols = packmol_system.nmols
    nfree = 0
    n_free_atoms = 0
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            nfree += st.number_of_molecules
            n_free_atoms += st.number_of_molecules * st.natoms
        end
    end

    # Build free_atom_mol (stable: depends only on system structure, never on positions)
    free_atom_mol = Vector{Int}(undef, n_free_atoms)
    ifree = 0
    imol_offset = 0
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            for i in 1:st.number_of_molecules
                imol = imol_offset + i
                for _ in 1:st.natoms
                    ifree += 1
                    free_atom_mol[ifree] = imol
                end
            end
        end
        imol_offset += st.number_of_molecules
    end

    x_vec = Vector{T}(undef, nfree * 2 * D)
    vaux = SPGBox.VAux(x_vec, zero(T))
    return MemoryBuffers{D,T,typeof(vaux)}(
        Vector{SVector{D,T}}(undef, natoms),
        zeros(T, nmols),
        Vector{SVector{D,T}}(undef, n_free_atoms),
        free_atom_mol,
        x_vec,
        vaux,
    )
end
