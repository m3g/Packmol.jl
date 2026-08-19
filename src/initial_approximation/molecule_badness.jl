#
# Compute per-molecule "badness" = constraint penalty + overlap with fixed molecules.
# Only considers constraint violations and distance violations against fixed atoms.
# Returns a vector of badness values indexed by molecule index.
#
# Two-phase computation:
#   Phase 1 (parallel): constraint penalties per molecule.
#   Phase 2 (sequential): fixed-atom overlaps via ParticleSystem2. pairwise! is
#     called once per molecule, updating only the x (free) cell list each time
#     while the y (fixed) cell list remains intact.
#
function molecule_badness(
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    fixed_sys::Union{Nothing,CellListMap.AbstractParticleSystem},
    tol::T;
    buffers::Union{Nothing,MemoryBuffers} = nothing,
) where {D,T}
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    nmols = packmol_system.nmols
    if isnothing(buffers)
        badness = zeros(T, nmols)
    else
        badness = buffers.badness
        fill!(badness, zero(T))
    end
    ntasks = Threads.nthreads()

    # Phase 1: constraint penalties (parallel over molecule ranges)
    _molecule_badness_constraints!(
        badness, packmol_system, atom_positions, Val(has_pbc), unitcell, unitcell_center, ntasks,
    )

    # Phase 2: fixed-atom overlaps — single pairwise! call over all free atoms.
    # xpositions is updated once with all free-atom positions; the fixed-atom
    # (ypositions) cell list is never rebuilt.
    if !isnothing(fixed_sys)
        if isnothing(buffers)
            n_free_atoms = sum(
                st.fixed.fixed ? 0 : st.number_of_molecules * st.natoms
                for st in packmol_system.structure_types
            )
            free_atoms = Vector{SVector{D,T}}(undef, n_free_atoms)
            free_atom_mol = Vector{Int}(undef, n_free_atoms)
            ifree = 0
            imol_offset = 0
            atom_offset = 0
            for st in packmol_system.structure_types
                if !st.fixed.fixed
                    for i in 1:st.number_of_molecules
                        imol = imol_offset + i
                        iat_first = atom_offset + (i - 1) * st.natoms + 1
                        for j in 0:st.natoms-1
                            ifree += 1
                            free_atoms[ifree] = atom_positions[iat_first + j]
                            free_atom_mol[ifree] = imol
                        end
                    end
                end
                imol_offset += st.number_of_molecules
                atom_offset += st.number_of_molecules * st.natoms
            end
        else
            free_atoms = buffers.free_atoms
            free_atom_mol = buffers.free_atom_mol
            _fill_free_atoms!(free_atoms, atom_positions, packmol_system)
        end
        # CellListMap's own cutoff may be larger than `tol` (grown to cap the
        # cell count for large/sparse boxes — see adjust_constraints.jl), so
        # pairs it finds aren't necessarily real violations; filter explicitly
        # rather than assuming every found pair has pair.d < tol.
        CellListMap.update!(fixed_sys; xpositions=free_atoms)
        pairwise!(
            (pair, out) -> begin
                pair.d < tol && (out.molecule_badness.value[free_atom_mol[pair.i]] += (pair.d - tol)^2)
                out
            end,
            fixed_sys,
        )
        badness .+= fixed_sys.output.molecule_badness.value
    end

    return badness
end

# HASPBC is a type parameter (via Val), not a runtime Bool, so the non-PBC
# specialization never compiles the wrap_to_center branch at all — see the
# comment on _constraint_fg! in interatomic_distance_fg.jl for why that matters.
function _molecule_badness_constraints!(
    badness::Vector{T},
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
    ntasks::Int,
) where {D,T,HASPBC}
    imol_offset = 0
    atom_offset = 0
    @sync for st in packmol_system.structure_types
        st_mol_offset = imol_offset
        st_atom_offset = atom_offset
        for (_, irange) in enumerate(chunks(1:st.number_of_molecules; n=ntasks))
            @spawn begin
                for i in irange
                    imol = st_mol_offset + i
                    iat_first = st_atom_offset + (i - 1) * st.natoms + 1
                    for j in 1:st.natoms
                        iat = iat_first + j - 1
                        x = atom_positions[iat]
                        if HASPBC
                            x = wrap_to_center(x, unitcell, unitcell_center)
                        end
                        for ic in packmol_system.atoms[iat].constraints
                            c = st.constraints[ic]
                            badness[imol] += constraint_penalty(c, x)
                        end
                    end
                end
            end
        end
        imol_offset += st.number_of_molecules
        atom_offset += st.number_of_molecules * st.natoms
    end
    return badness
end
