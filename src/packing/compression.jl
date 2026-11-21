#
# Packing by compression (`use_compression` keyword): the system is packed at a lower density,
# with its constraints and periodic cell expanded, and then compressed back in steps,
# each packed from the previous one. Packings too dense to converge from a random
# placement (long rigid molecules, or polymers, that jam) converge this way.
#
# At each step the geometry is scaled by a factor `s < 1` about the center of the cell
# (the origin without periodic boundary conditions): the constraints, the cell, and the
# positions of the free molecules, each moved rigidly (its center of mass scaled; a flexible
# molecule by the center of its fragments, such that the links between them are kept), and
# the fixed molecules, scaled with their radii (see `_FixedScaling`). The step
# (`1 - s`) starts at the `initial step` given, and grows by 1.5 after two steps that
# converge, up to the initial step. When a step does not converge, the system "breathes":
# the packing of the failed step (with the rearrangements it made) is expanded to one step
# below the density of the last converged one, and packed there; the compression resumes
# from it, with half the step. If the relaxed packing does not converge either, or after
# `_compress_max_breaths` breaths, the compression goes back to the last converged packing
# with half the step. Going back alone stopped at ~90% of the density in dense systems
# (lignin particle at 1.3 g/cm³, 2 of 3 runs): the failed step was discarded, and the
# next steps failed at the same configuration.
#

# The loops of each step of the compression
const _compress_step_loops = 100
# The smallest step: below it, the compression stops
const _compress_min_step = 5e-4
# The largest number of breaths (relaxations of a failed step) of a compression
const _compress_max_breaths = 20
# The largest number of further expansions when the packing at low density does not converge
const _compress_max_expansions = 3

#
# The constraints scaled by `s` about the point `o`: the region of `c` mapped by
# `x -> o + s (x - o)`
#
_scaled_constraint(c::Sphere{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Sphere{P,T}(o + s * (c.center - o), s * c.radius, c.weight)
_scaled_constraint(c::Box{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Box{P,T}(o + s * (c.center - o), s * c.sides, c.weight)
_scaled_constraint(c::Cube{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Cube{P,T}(o + s * (c.center - o), s * c.side, c.weight)
_scaled_constraint(c::Cylinder{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Cylinder{P,T}(o + s * (c.center - o), c.axis, s * c.radius, s * c.length, c.weight)
_scaled_constraint(c::Ellipsoid{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Ellipsoid{P,T}(o + s * (c.center - o), s * c.a, s * c.b, s * c.c, c.scale, c.weight)
# `normal ⋅ x = d`, with `x = o + (x′ - o) / s`
_scaled_constraint(c::Plane{P,T}, s::T, o::SVector{3,T}) where {P,T} =
    Plane{P,T}(c.normal, s * c.d + (1 - s) * dot(c.normal, o), c.weight)
_scaled_constraint(c, ::Any, ::Any) = throw(ArgumentError("""\n
    use_compression: the constraints of type $(nameof(typeof(c))) cannot be scaled (supported: boxes,
    cubes, spheres, cylinders, ellipsoids and planes).

"""))

# Whether the geometry of the system can be scaled (an error otherwise)
function _check_compress(sys::PackmolSystem{D,T}) where {D,T}
    D == 3 || throw(ArgumentError("""\n
        use_compression is supported in three dimensions only.

    """))
    isnothing(sys.voxel_place) || throw(ArgumentError("""\n
        use_compression is not supported with voxel_place.

    """))
    isnothing(sys.restart_from) || throw(ArgumentError("""\n
        use_compression is not supported with restart_from.

    """))
    o = zero(SVector{D,T})
    for st in sys.structure_types, c in st.constraints
        _constraint_dispatch(_scaled_constraint, c, one(T), o)
    end
    return nothing
end

# The center of the scaling: the center of the cell, or the origin
_compress_center(sys::PackmolSystem{D,T}) where {D,T} =
    isnothing(sys.unitcell_center) ? zero(SVector{D,T}) : SVector{D,T}(sys.unitcell_center)

#
# Scales the geometry of the system by `s` about `o`: the constraints, the cell, and (if
# `molecules`) the positions of the free molecules, each moved rigidly. The fixed
# molecules are scaled apart (`_scale_fixed!`), from their original coordinates
#
function _scale_system!(sys::PackmolSystem{D,T}, s::T, o::SVector{D,T}; molecules::Bool=true) where {D,T}
    # (the polymers as rigid molecules have constraints of their own, for their placement)
    for st in Iterators.flatten((sys.structure_types, (rst for p in sys.polymers for rst in p.rigid_types)))
        for i in eachindex(st.constraints)
            st.constraints[i] = _constraint_dispatch(_scaled_constraint, st.constraints[i], s, o)
        end
    end
    unitcell = sys.unitcell
    if !isnothing(unitcell)
        sys.unitcell = s .* unitcell
    end
    molecules || return sys
    topo = sys.polymer_topology
    imol = 0
    for st in sys.structure_types, _ in 1:st.number_of_molecules
        imol += 1
        (st.fixed.fixed || _is_monomer(topo, imol)) && continue
        mp = sys.molecule_positions[imol]
        sys.molecule_positions[imol] = MoleculePosition(o + s * (mp.cm - o), mp.angles)
    end
    # Each polymer moved by the shift of its center (its chain made whole first: the
    # centers of mass of the monomers are wrapped into the cell one by one)
    mol_offset, _ = _type_offsets(sys)
    has_pbc, unitcell_matrix, _ = _typed_unitcell(sys)
    for p in sys.polymers
        for (c, j) in _polymer_copies(p)
            mols = [_monomer_molecule(p, mol_offset, c, m, j) for m in 1:_nmonomers(p)]
            previous = sys.molecule_positions[mols[1]].cm
            center = zero(SVector{D,T})
            for imol in mols
                cm = sys.molecule_positions[imol].cm
                x = has_pbc ? previous + _min_image(cm - previous, unitcell_matrix) : cm
                center += x
                previous = x
            end
            center /= length(mols)
            shift = (s - 1) * (center - o)
            for imol in mols
                mp = sys.molecule_positions[imol]
                sys.molecule_positions[imol] = MoleculePosition(mp.cm + shift, mp.angles)
            end
        end
    end
    return sys
end

#
# The fixed molecules are scaled with the geometry: the coordinates of their atoms about
# the center of the scaling, and the radii of their atoms, such that the constraints placed
# relative to them keep their relation to them (the membrane of the Zika virion, packed
# by compression, crossed by the transmembrane segments of its fixed proteins, could not
# be compressed with the proteins at their place: it stopped at a third of the density).
# The original coordinates of each fixed type and radii of the fixed atoms, from which
# the fixed molecules are set at each scale (restored exactly at the scale of the input).
#
struct _FixedScaling{D,T}
    types::Vector{Int}
    coordinates::Vector{Vector{SVector{D,T}}}
    radii::Vector{Vector{T}}
    atoms::Vector{Int}
    atom_radii::Vector{T}
end

function _FixedScaling(sys::PackmolSystem{D,T}) where {D,T}
    types = Int[]
    coordinates = Vector{SVector{D,T}}[]
    radii = Vector{T}[]
    for (t, st) in enumerate(sys.structure_types)
        st.fixed.fixed || continue
        push!(types, t)
        push!(coordinates, copy(st.reference_coordinates))
        push!(radii, copy(st.radii))
    end
    atoms = Int[]
    atom_radii = T[]
    for (i, a) in enumerate(sys.atoms)
        sys.structure_types[a.structure_type_index].fixed.fixed || continue
        push!(atoms, i)
        push!(atom_radii, a.radius)
    end
    return _FixedScaling{D,T}(types, coordinates, radii, atoms, atom_radii)
end

# The fixed molecules at the scale `k` of the geometry, about `o`
function _scale_fixed!(sys::PackmolSystem{D,T}, state::_FixedScaling{D,T}, k::T, o::SVector{D,T}) where {D,T}
    exact = abs(k - one(T)) <= sqrt(eps(T))
    for (n, t) in enumerate(state.types)
        st = sys.structure_types[t]
        for i in eachindex(st.reference_coordinates)
            x0 = state.coordinates[n][i]
            st.reference_coordinates[i] = exact ? x0 : o + k * (x0 - o)
        end
        for i in eachindex(st.radii)
            st.radii[i] = exact ? state.radii[n][i] : k * state.radii[n][i]
        end
    end
    for (n, i) in enumerate(state.atoms)
        a = sys.atoms[i]
        r0 = state.atom_radii[n]
        sys.atoms[i] = AtomData(a.molecule_index, a.structure_type_index, exact ? r0 : k * r0, a.constraints)
    end
    return sys
end

#
# The packing by compression: the system packed with its geometry expanded by
# `fraction^(-1/3)` (`fraction` of the density), then compressed in steps (see above).
# `run_options` are the keyword arguments of `_packmol`.
#
function _packmol_compressed(sys::PackmolSystem{D,T}; run_options...) where {D,T}
    fraction, initial_step = something(sys.use_compression)
    _check_compress(sys)
    o = _compress_center(sys)
    expansion = T(fraction)^(-one(T) / 3)
    _println()
    @_printf("  Packing at %.1f%% of the density (the geometry expanded by %.4f), then compressing\n",
        100 * fraction, expansion)
    # The constraints and the cell expanded (the molecules are placed by the packing)
    _scale_system!(sys, expansion, o; molecules=false)
    fixed_state = _FixedScaling(sys)
    _scale_fixed!(sys, fixed_state, expansion, o)
    # The initial placement, at the low density
    sys.check_initial_placement && return _packmol(sys; run_options...)
    output_file = sys.output_file
    sys.output_file = ""
    scale = expansion  # the current scale of the geometry, relative to that of the input
    try
        _packmol(sys; run_options...)
        # If the packing at low density does not converge, the geometry is expanded
        # further (half the density each time, the molecules moved with it) and the
        # system packed again from there, up to `_compress_max_expansions` times: the
        # low-density packing of the polymer-capped gold nanoparticle converged in some
        # runs and not in others
        nexpansions = 0
        while sys.status != :packing_ready && nexpansions < _compress_max_expansions
            nexpansions += 1
            factor = T(2)^(one(T) / 3)
            _scale_system!(sys, factor, o)
            scale *= factor
            _scale_fixed!(sys, fixed_state, scale, o)
            _println()
            @_printf("  Not converged: packing again at %.1f%% of the density (the geometry expanded by %.4f)\n",
                100 / scale^3, scale)
            _packmol(sys; restart=true, run_options...)
        end
        if sys.status != :packing_ready
            _println("  The packing at low density did not converge: the compression was not done.")
            return sys
        end
        # The compression steps, without the report of each packing
        screen_log = _screen_log[]
        _screen_log[] = false
        step = T(initial_step)
        nconverged = 0
        nsteps = 0
        seed = sys.seed
        tstart = time()
        nbreaths = 0
        try
            while scale > one(T) + sqrt(eps(T))
                s = max(one(T) - step, one(T) / scale)
                saved_positions = copy(sys.molecule_positions)
                saved_constraints = [copy(st.constraints) for st in sys.structure_types]
                saved_unitcell = sys.unitcell
                _scale_system!(sys, s, o)
                _scale_fixed!(sys, fixed_state, scale * s, o)
                nsteps += 1
                seed += 1
                t0 = time()
                _packmol(sys; restart=true, nloop=_compress_step_loops, seed, run_options...)
                converged = sys.status == :packing_ready
                if screen_log
                    _screen_log[] = true
                    @_printf("  Compression step %4d: density %7.2f%% of the target, %s (%.1f s)\n",
                        nsteps, 100 / (scale * s)^3, converged ? "converged" : "not converged", time() - t0)
                    _screen_log[] = false
                end
                if converged
                    scale *= s
                    nconverged += 1
                    if nconverged >= 2
                        step = min(T(initial_step), T(1.5) * step)
                        nconverged = 0
                    end
                elseif nbreaths < _compress_max_breaths
                    # the failed packing relaxed: expanded to one step below the density
                    # of the last converged one (keeping the rearrangements of the failed
                    # step), and packed; the compression resumes from it with half the
                    # step (or from the last converged packing, if it does not converge)
                    nbreaths += 1
                    relaxed_scale = scale / (one(T) - step)
                    _scale_system!(sys, relaxed_scale / (scale * s), o)
                    _scale_fixed!(sys, fixed_state, relaxed_scale, o)
                    t0 = time()
                    _packmol(sys; restart=true, nloop=_compress_step_loops, seed, run_options...)
                    relaxed = sys.status == :packing_ready
                    if screen_log
                        _screen_log[] = true
                        @_printf("  Relaxation %4d:        density %7.2f%% of the target, %s (%.1f s)\n",
                            nbreaths, 100 / relaxed_scale^3, relaxed ? "converged" : "not converged", time() - t0)
                        _screen_log[] = false
                    end
                    if relaxed
                        scale = relaxed_scale
                    else
                        copyto!(sys.molecule_positions, saved_positions)
                        for (st, cs) in zip(sys.structure_types, saved_constraints)
                            copyto!(st.constraints, cs)
                        end
                        sys.unitcell = saved_unitcell
                        _scale_fixed!(sys, fixed_state, scale, o)
                    end
                    step /= 2
                    nconverged = 0
                    step < _compress_min_step && break
                else
                    # back to the last converged packing, with half the step
                    copyto!(sys.molecule_positions, saved_positions)
                    for (st, cs) in zip(sys.structure_types, saved_constraints)
                        copyto!(st.constraints, cs)
                    end
                    sys.unitcell = saved_unitcell
                    _scale_fixed!(sys, fixed_state, scale, o)
                    step /= 2
                    nconverged = 0
                    step < _compress_min_step && break
                end
            end
        finally
            _screen_log[] = screen_log
            _scale_fixed!(sys, fixed_state, scale, o)
        end
        reached = scale <= one(T) + sqrt(eps(T))
        sys.status = reached ? :packing_ready : :failed_packing
        _println()
        if reached
            @_printf("  Compressed to the density of the input in %d steps (%.1f s).\n", nsteps, time() - tstart)
        else
            @_printf("  WARNING: the compression stopped at %.2f%% of the density of the input (a step\n", 100 / scale^3)
            _println("  smaller than the minimum did not converge). The geometry of the output is that of")
            _println("  this density.")
        end
    finally
        sys.output_file = output_file
    end
    if !isempty(sys.output_file)
        write_output(sys)
        _println("  Solution written to file: ", sys.output_file)
    end
    return sys
end
