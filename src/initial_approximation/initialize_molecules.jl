#
# Initialize molecule positions randomly within constraint bounding box
# (or within the unit cell for PBC), and center reference coordinates
# at the origin (required for chain rule).
#
# This is just the first rough placement — no constraint checking is done here.
# The subsequent steps will estimate bounds from the best molecules and
# re-initialize within those bounds.
#

#
# For the common case of a structure type constrained entirely by whole-molecule
# Box/Cube/Sphere constraints (i.e. every atom references the same full set of
# constraints — no per-atom sub-constraints), derive a direct CM sampling region
# from the geometry of the "Inside" constraints (their bounding-box intersection,
# shrunk by the molecule's maximum reference-coordinate extent so a randomly
# rotated molecule tends to stay inside), and reject/retry samples against the
# full constraint set (including any "Outside" constraints, e.g. a spherical
# shell built from `inside sphere` + `outside sphere`).
#
# This matches the original Fortran Packmol, which samples molecule centers
# directly and roughly uniformly within the true constraint region. Without it,
# `initialize_molecules!` falls back to a huge generic bootstrap box (`sidemax`)
# followed by unconstrained gradient-only constraint fitting with no inter-
# molecule repulsion term — which pulls every molecule straight to the *nearest*
# feasible point and stops as soon as it crosses the boundary (the penalty and
# its gradient vanish once inside), so molecules starting outside the region all
# pile up right at its surface instead of spreading through the volume — e.g.
# for a spherical shell, everything collapses onto the outer sphere.
#
# Returns `nothing` if no direct region can be derived (per-atom sub-constraints,
# no Box/Cube/Sphere Inside constraint to bound the region, or a shrink margin
# too large relative to the region), in which case the caller falls back to the
# sidemax bootstrap.
#
function _inside_bound(c::Union{Box{Inside},Cube{Inside},Sphere{Inside}}, margin::T) where {T}
    center, sides = if c isa Box{Inside}
        c.center, c.sides
    elseif c isa Cube{Inside}
        c.center, SVector{3,T}(c.side, c.side, c.side)
    else
        c.center, SVector{3,T}(2c.radius, 2c.radius, 2c.radius)
    end
    half = sides ./ 2
    lo, hi = center .- half, center .+ half
    m = SVector{3,T}(margin, margin, margin)
    lo_shrunk, hi_shrunk = lo .+ m, hi .- m
    # If the margin is too large relative to the region, the shrink would
    # invert the range; fall back to the raw (unshrunk) region instead.
    return all(lo_shrunk .<= hi_shrunk) ? (lo_shrunk, hi_shrunk) : (lo, hi)
end
_inside_bound(::Constraint, margin) = nothing

function _direct_cm_region(st::StructureType{D,T}, max_extent::T) where {D,T}
    n = length(st.constraints)
    n == 0 && return nothing
    all(ac -> ac == 1:n, st.atom_constraints) || return nothing
    lo = SVector{3,T}(ntuple(_ -> T(-Inf), 3))
    hi = SVector{3,T}(ntuple(_ -> T(Inf), 3))
    has_inside = false
    for c in st.constraints
        bounds = _inside_bound(c, max_extent)
        bounds === nothing && continue
        clo, chi = bounds
        lo, hi = max.(lo, clo), min.(hi, chi)
        has_inside = true
    end
    (has_inside && all(lo .<= hi)) || return nothing
    return lo, hi
end

#
# Margin-adjusted copies of a constraint, used only for the CM-point rejection
# check during initial placement: "Inside" regions are shrunk and "Outside"
# regions are grown by the molecule's max reference-coordinate extent, so that
# a CM satisfying the adjusted constraint leaves room for the whole (randomly
# rotated) molecule to stay within the real constraint. Falls back to the
# unshrunk constraint if the margin would invert an "Inside" region.
#
_margin_adjust(c::Sphere{Inside,T}, margin::T) where {T} =
    (r = c.radius - margin; Sphere{Inside,T}(center=c.center, radius=(r > 0 ? r : c.radius), weight=c.weight))
_margin_adjust(c::Sphere{Outside,T}, margin::T) where {T} =
    Sphere{Outside,T}(center=c.center, radius=c.radius + margin, weight=c.weight)
_margin_adjust(c::Box{Inside,T}, margin::T) where {T} =
    (s = c.sides .- 2margin; Box{Inside,T}(center=c.center, sides=(all(s .> 0) ? s : c.sides), weight=c.weight))
_margin_adjust(c::Box{Outside,T}, margin::T) where {T} =
    Box{Outside,T}(center=c.center, sides=c.sides .+ 2margin, weight=c.weight)
_margin_adjust(c::Cube{Inside,T}, margin::T) where {T} =
    (s = c.side - 2margin; Cube{Inside,T}(center=c.center, side=(s > 0 ? s : c.side), weight=c.weight))
_margin_adjust(c::Cube{Outside,T}, margin::T) where {T} =
    Cube{Outside,T}(center=c.center, side=c.side + 2margin, weight=c.weight)
_margin_adjust(c::Constraint, margin) = c

# Best-of-`max_tries` uniform sample within the AABB [lo, hi], scored by total
# penalty against `constraints`; returns as soon as a zero-penalty sample is found.
function _sample_in_region(rng, lo::SVector{3,T}, hi::SVector{3,T}, constraints, max_tries::Int) where {T}
    extent = hi .- lo
    draw() = SVector{3,T}(ntuple(d -> lo[d] + rand(rng, T) * extent[d], 3))
    best = draw()
    best_pen = sum(c -> constraint_penalty(c, best), constraints; init=zero(T))
    tries = 1
    while best_pen > zero(T) && tries < max_tries
        cm = draw()
        pen = sum(c -> constraint_penalty(c, cm), constraints; init=zero(T))
        if pen < best_pen
            best_pen, best = pen, cm
        end
        tries += 1
    end
    return best
end

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
    # Packmol's sidemax approach) as a fallback for structure types whose
    # constraint geometry doesn't admit a direct CM sampling region (see
    # _direct_cm_region above). The constraint-only optimization that follows
    # will move molecules to feasible positions regardless of constraint type.
    sidemax = T(DEFAULT_SIDEMAX)
    max_region_tries = 30

    # Compute molecule index offset for each structure type so threads
    # can determine the correct slot without a shared counter.
    imol_offset = 0
    @sync for st in packmol_system.structure_types
        st_offset = imol_offset
        direct_region, direct_constraints = if !has_pbc && !st.fixed.fixed
            max_extent = isempty(st.reference_coordinates) ? zero(T) : maximum(norm, st.reference_coordinates)
            region = _direct_cm_region(st, max_extent)
            if region === nothing
                nothing, nothing
            else
                region, [_margin_adjust(c, max_extent) for c in st.constraints]
            end
        else
            nothing, nothing
        end
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
                        elseif direct_region !== nothing
                            lo, hi = direct_region
                            cm3 = _sample_in_region(task_rng, lo, hi, direct_constraints, max_region_tries)
                            cm = SVector{D,T}(ntuple(d -> cm3[d], D))
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
