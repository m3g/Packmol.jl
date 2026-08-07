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
# For the common case of a structure type constrained by a single
# whole-molecule Inside box/cube (i.e. every atom references exactly that
# one constraint — no per-atom sub-constraints), derive a direct CM
# sampling region from the constraint geometry itself, shrunk by the
# molecule's maximum reference-coordinate extent so a randomly rotated
# molecule tends to stay inside.
#
# This matches the original Fortran Packmol, which samples molecule
# centers directly and uniformly within the true constraint region for
# this case. Without it, `initialize_molecules!` falls back to a huge
# generic bootstrap box (`sidemax`) followed by unconstrained gradient-only
# constraint fitting with no inter-molecule repulsion term — which herds
# many molecules toward similar feasible regions/box faces, producing
# initial overlap that grows worse than linearly with molecule count.
#
# Returns `nothing` if no such direct region can be derived (multiple/
# per-atom constraints, non-box/cube constraint, or Outside placement),
# in which case the caller should fall back to the sidemax bootstrap.
#
function _direct_cm_bounds(st::StructureType{D,T}, max_extent::T) where {D,T}
    (length(st.constraints) == 1 && all(==([1]), st.atom_constraints)) || return nothing
    c = st.constraints[1]
    lo, hi = if c isa Box{Inside}
        c.center .- c.sides ./ 2, c.center .+ c.sides ./ 2
    elseif c isa Cube{Inside}
        c.center .- c.side / 2, c.center .+ c.side / 2
    else
        return nothing
    end
    margin = SVector{D,T}(ntuple(_ -> max_extent, D))
    lo_shrunk, hi_shrunk = lo .+ margin, hi .- margin
    # If the molecule is larger than the box, the margin would invert the
    # range; fall back to sampling within the raw (unshrunk) box instead
    # of an empty region.
    return all(lo_shrunk .<= hi_shrunk) ? (lo_shrunk, hi_shrunk) : (lo, hi)
end

#
# Analogous direct-sampling shortcut for a whole-molecule Inside sphere
# constraint: sample the CM uniformly within the ball (shrunk by the
# molecule's max reference-coordinate extent), instead of falling back to
# the sidemax bootstrap box. For a sphere, the constraint-only gradient
# fitting used by the bootstrap path pulls every molecule straight to the
# nearest feasible point and stops as soon as it crosses the boundary
# (the penalty and its gradient vanish once inside) — so molecules starting
# outside the sphere all pile up just inside its surface instead of
# spreading through the volume.
#
function _direct_cm_sphere(st::StructureType{D,T}, max_extent::T) where {D,T}
    (length(st.constraints) == 1 && all(==([1]), st.atom_constraints)) || return nothing
    c = st.constraints[1]
    c isa Sphere{Inside} || return nothing
    radius_shrunk = c.radius - max_extent
    radius = radius_shrunk > 0 ? radius_shrunk : c.radius
    return c.center, radius
end

# Uniformly sample a point within the D-dimensional ball of given center/radius.
function _random_point_in_ball(rng, center::SVector{D,T}, radius::T) where {D,T}
    dir = SVector{D,T}(ntuple(_ -> randn(rng, T), D))
    n = norm(dir)
    dir = n > 0 ? dir / n : SVector{D,T}(ntuple(d -> d == 1 ? one(T) : zero(T), D))
    r = radius * rand(rng, T)^(one(T) / D)
    return center + r * dir
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
    # _direct_cm_bounds above). The constraint-only optimization that follows
    # will move molecules to feasible positions regardless of constraint type.
    sidemax = T(DEFAULT_SIDEMAX)

    # Compute molecule index offset for each structure type so threads
    # can determine the correct slot without a shared counter.
    imol_offset = 0
    @sync for st in packmol_system.structure_types
        st_offset = imol_offset
        direct_bounds, direct_sphere = if !has_pbc && !st.fixed.fixed
            max_extent = isempty(st.reference_coordinates) ? zero(T) : maximum(norm, st.reference_coordinates)
            bounds = _direct_cm_bounds(st, max_extent)
            bounds, bounds === nothing ? _direct_cm_sphere(st, max_extent) : nothing
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
                        elseif direct_bounds !== nothing
                            lo, hi = direct_bounds
                            cm = SVector{D,T}(ntuple(d -> lo[d] + rand(task_rng, T) * (hi[d] - lo[d]), D))
                        elseif direct_sphere !== nothing
                            sphere_center, sphere_radius = direct_sphere
                            cm = _random_point_in_ball(task_rng, sphere_center, sphere_radius)
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
