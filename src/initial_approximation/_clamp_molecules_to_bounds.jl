#
# Clamp molecules whose CMs are far outside the per-type constraint bounds.
# After initial constraint adjustment, some molecules may still be at sidemax
# coordinates if adjustment did not converge. We move those molecules randomly
# inside the constraint bounds to prevent a huge bounding box.
#
function _clamp_molecules_to_bounds!(
    packmol_system::PackmolSystem{D,T},
    free_mol_indices::Vector{Int},
    cm_min::Vector{SVector{D,T}},
    cm_max::Vector{SVector{D,T}},
    RNG,
) where {D,T}
    imol = 0
    nclamped = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        lo = cm_min[ist]
        hi = cm_max[ist]
        extent = hi - lo
        has_valid_bounds = all(extent .> zero(T))
        for _ in 1:st.number_of_molecules
            imol += 1
            st.fixed.fixed && continue
            mp = packmol_system.molecule_positions[imol]
            cm = mp.cm
            if has_valid_bounds && any(cm .< lo) || any(cm .> hi)
                # Move to a random position within the bounds
                new_cm = lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
                packmol_system.molecule_positions[imol] = MoleculePosition(new_cm, mp.angles)
                nclamped += 1
            end
        end
    end
    if nclamped > 0
        println("  Clamped $nclamped molecules to constraint bounds.")
    end
    return nothing
end