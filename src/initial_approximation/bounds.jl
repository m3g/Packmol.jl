#
# Bounds of the positions and of the centers of mass of the molecules.
#

#
# Bounding box computed from actual atom positions.
# Used after molecule placement to determine the region for CellListMap.
#
function compute_bounding_box(atom_positions::Vector{SVector{D,T}}) where {D,T}
    lo = SVector{D,T}(ntuple(_ -> typemax(T), D))
    hi = SVector{D,T}(ntuple(_ -> typemin(T), D))
    for x in atom_positions
        lo = min.(lo, x)
        hi = max.(hi, x)
    end
    return lo, hi
end

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
                if precision == typemax(T) || constraint_penalty_sum(mp, st, packmol_system) ≤ precision
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
            half_extent = _cm_half_extent(packmol_system)
            lo = max.(lo, min.(cm_min[ist], packmol_system.unitcell_center .- half_extent))
            hi = min.(hi, max.(cm_max[ist], packmol_system.unitcell_center .+ half_extent))
        end
        cm_min[ist] = lo
        cm_max[ist] = hi
    end
    return cm_min, cm_max
end

#
# The bounds of the variables of `n` molecules in the constraint fits of the initial
# approximation and of movebad: under periodic boundary conditions, the centers of mass
# in the box of the cell (`_cm_half_extent`, as in the packing loop), so that they don't
# drift along the directions that no constraint bounds (a molecule above a plane, along
# the plane); without PBC, and for the angles, none (infinite bounds, for the calls to
# have a single type).
#
# They replace the implicit confinement of every atom to the cell, which pushed the
# molecules longer than the tolerance away from the faces of the cell in the fits (the
# lipids of a bilayer crowded in the middle of the cell, 15 of 60 lipids of a leaflet in
# the outer quarter of the square, against 26 for a uniform distribution).
#
function _fit_bounds(packmol_system::PackmolSystem{D,T}, n::Int) where {D,T}
    inf = SVector{D,T}(ntuple(_ -> T(Inf), D))
    lo, hi = if isnothing(packmol_system.unitcell)
        -inf, inf
    else
        h = _cm_half_extent(packmol_system)
        c = packmol_system.unitcell_center::SVector{D,T}
        c - h, c + h
    end
    lower = collect(reinterpret(T, fill(MoleculePosition(lo, -inf), n)))
    upper = collect(reinterpret(T, fill(MoleculePosition(hi, inf), n)))
    return lower, upper
end

#
# Clamp molecules whose CMs are far outside the per-type constraint bounds.
# After initial constraint adjustment, some molecules may still be at sidemax
# coordinates if adjustment did not converge. We move those molecules randomly
# inside the constraint bounds to prevent a huge bounding box. Molecules that
# satisfy their constraints are kept where they are: the bounds are a box around
# the region of the constraints, and a random point in it may not satisfy them.
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
            if has_valid_bounds && (any(cm .< lo) || any(cm .> hi)) &&
               constraint_penalty_sum(mp, st, packmol_system) > packmol_system.constraint_precision
                # Move to a random position within the bounds
                new_cm = lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
                packmol_system.molecule_positions[imol] = MoleculePosition(new_cm, mp.angles)
                nclamped += 1
            end
        end
    end
    if nclamped > 0
        _println("  Clamped $nclamped molecules to constraint bounds.")
    end
    return nothing
end
