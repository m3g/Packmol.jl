#
# Compute per-structure-type center-of-mass bounds from current molecule positions.
# Returns (cm_min, cm_max) vectors indexed by structure type.
# Used after constraint fitting to determine where molecules should be placed.
#
function compute_cm_bounds(packmol_system::PackmolSystem{D,T}; precision::T=typemax(T)) where {D,T}
    ntypes = length(packmol_system.structure_types)
    cm_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    cm_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]
    # Track per-type fallback bounds (all molecules), used when no molecule
    # satisfies the constraint precision for a given structure type.
    fb_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    fb_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]
    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        for _ in 1:st.number_of_molecules
            imol += 1
            if !st.fixed.fixed
                mp = packmol_system.molecule_positions[imol]
                cm = mp.cm
                fb_min[ist] = min.(fb_min[ist], cm)
                fb_max[ist] = max.(fb_max[ist], cm)
                if precision == typemax(T) || constraint_penalty_sum(mp, st) ≤ precision
                    cm_min[ist] = min.(cm_min[ist], cm)
                    cm_max[ist] = max.(cm_max[ist], cm)
                end
            end
        end
    end
    # Fallback: for any type where no molecule satisfied the precision filter,
    # use the unconstrained bounds (all molecules).
    for ist in 1:ntypes
        if any(cm_min[ist] .== typemax(T))
            cm_min[ist] = fb_min[ist]
            cm_max[ist] = fb_max[ist]
        end
    end
    return cm_min, cm_max
end