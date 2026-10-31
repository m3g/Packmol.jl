#
# Compute per-structure-type center-of-mass bounds from current molecule positions.
# Returns (cm_min, cm_max) vectors indexed by structure type.
# Used after constraint fitting to determine where molecules should be placed.
#
# With `extrapolate = true` (the default), the range of each type is widened by
# range/(n - 1) on each side, n being the number of molecules it was computed
# from: the observed range of n centers underestimates the region they can
# occupy (for n uniform points in an interval, the expected gap at each end is
# 1/(n + 1) of it). Without this, a type of few molecules (e.g. one of several
# conformations of a lipid, 6 copies each) is placed only within the box
# spanned by its own few centers, a random part of its region, and the gaps
# elsewhere are never filled. Under PBC, the widened bounds are clamped to the
# cell.
#
function compute_cm_bounds(
    packmol_system::PackmolSystem{D,T};
    precision::T=typemax(T),
    extrapolate::Bool=true,
) where {D,T}
    ntypes = length(packmol_system.structure_types)
    cm_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    cm_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]
    # Track per-type fallback bounds (all molecules), used when no molecule
    # satisfies the constraint precision for a given structure type.
    fb_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    fb_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]
    n_in = zeros(Int, ntypes)
    n_fb = zeros(Int, ntypes)
    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        for _ in 1:st.number_of_molecules
            imol += 1
            if !st.fixed.fixed
                mp = packmol_system.molecule_positions[imol]
                cm = mp.cm
                fb_min[ist] = min.(fb_min[ist], cm)
                fb_max[ist] = max.(fb_max[ist], cm)
                n_fb[ist] += 1
                if precision == typemax(T) || constraint_penalty_sum(mp, st) ≤ precision
                    cm_min[ist] = min.(cm_min[ist], cm)
                    cm_max[ist] = max.(cm_max[ist], cm)
                    n_in[ist] += 1
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
            n_in[ist] = n_fb[ist]
        end
    end
    extrapolate || return cm_min, cm_max
    unitcell = packmol_system.unitcell
    for ist in 1:ntypes
        n = n_in[ist]
        n < 2 && continue
        margin = (cm_max[ist] .- cm_min[ist]) ./ (n - 1)
        lo = cm_min[ist] .- margin
        hi = cm_max[ist] .+ margin
        if !isnothing(unitcell)
            # Axis-aligned half-extent of the cell, as the bounds of `packmol`
            half_extent = SVector{D,T}(ntuple(i -> sum(abs(unitcell[i, j]) for j in 1:D) / 2, D))
            lo = max.(lo, min.(cm_min[ist], packmol_system.unitcell_center .- half_extent))
            hi = min.(hi, max.(cm_max[ist], packmol_system.unitcell_center .+ half_extent))
        end
        cm_min[ist] = lo
        cm_max[ist] = hi
    end
    return cm_min, cm_max
end