#
# Initialize molecule positions randomly within constraint bounding box
# (or within the unit cell for PBC), and center reference coordinates
# at the origin (required for chain rule).
#
# This is just the first rough placement — no constraint checking is done here.
# The subsequent steps will estimate bounds from the best molecules and
# re-initialize within those bounds.
#
function initialize_molecules!(packmol_system::PackmolSystem{D,T}, RNG) where {D,T}
    # Center reference coordinates at origin (required for chain rule)
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            cm = mean(st.reference_coordinates)
            st.reference_coordinates .-= Ref(cm)
        end
    end

    # Determine the placement region
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
    end

    # Initial placement: use a large box centered at origin (following Fortran
    # Packmol's sidemax approach). The constraint-only optimization that follows
    # will move molecules to feasible positions regardless of constraint type.
    sidemax = T(DEFAULT_SIDEMAX)

    # Compute molecule index offset for each structure type so threads
    # can determine the correct slot without a shared counter.
    imol_offset = 0
    @sync for st in packmol_system.structure_types
        st_offset = imol_offset
        for (_, irange) in enumerate(chunks(1:st.number_of_molecules; n=Threads.nthreads()))
            task_seed = rand(RNG, UInt64)
            @spawn begin
                task_rng = typeof(RNG)(task_seed)
                for i in irange
                    imol_local = st_offset + i
                    if st.fixed.fixed
                        packmol_system.molecule_positions[imol_local] = st.fixed.position
                    else
                        if has_pbc
                            frac = SVector{D,T}(ntuple(_ -> rand(task_rng, T) - T(0.5), D))
                            cm = SVector{D,T}(uc * frac) + center
                        else
                            cm = SVector{D,T}(ntuple(_ -> sidemax * (T(2) * rand(task_rng, T) - one(T)), D))
                        end
                        angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(task_rng, T), D))
                        packmol_system.molecule_positions[imol_local] = MoleculePosition(cm, angles)
                    end
                end
            end
        end
        imol_offset += st.number_of_molecules
    end
    return packmol_system
end
