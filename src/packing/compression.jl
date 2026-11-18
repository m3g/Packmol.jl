#
# Packing by compression (`compress` keyword): the system is packed at a lower density,
# with its constraints and periodic cell expanded, and then compressed back in steps,
# each packed from the previous one. Packings too dense to converge from a random
# placement (long rigid molecules, or polymers, that jam) converge this way.
#
# At each step the geometry is scaled by a factor `s < 1` about the center of the cell
# (the origin without periodic boundary conditions): the constraints, the cell, and the
# positions of the molecules, each moved rigidly (its center of mass scaled; a polymer
# by the center of its monomers, such that the links between them are kept). The step
# (`1 - s`) starts at the `initial step` given, is halved when a step does not converge
# (the packing goes back to the last converged one), and grows by 1.5 after two steps
# that converge, up to the initial step.
#

# The loops of each step of the compression
const _compress_step_loops = 100
# The smallest step: below it, the compression stops
const _compress_min_step = 5e-4

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
    compress: the constraints of type $(nameof(typeof(c))) cannot be scaled (supported: boxes,
    cubes, spheres, cylinders, ellipsoids and planes).

"""))

# Whether the geometry of the system can be scaled (an error otherwise)
function _check_compress(sys::PackmolSystem{D,T}) where {D,T}
    D == 3 || throw(ArgumentError("""\n
        compress is supported in three dimensions only.

    """))
    any(st -> st.fixed.fixed, sys.structure_types) && throw(ArgumentError("""\n
        compress is not supported with fixed molecules: they would have to move with the
        compression. Pack the free molecules around them in a second step.

    """))
    isnothing(sys.restart_from) || throw(ArgumentError("""\n
        compress is not supported with restart_from.

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
# `molecules`) the positions of the free molecules, each moved rigidly
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
    for imol in 1:sys.nmols
        _is_monomer(topo, imol) && continue
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
# The packing by compression: the system packed with its geometry expanded by
# `fraction^(-1/3)` (`fraction` of the density), then compressed in steps (see above).
# `run_options` are the keyword arguments of `_packmol`.
#
function _packmol_compressed(sys::PackmolSystem{D,T}; run_options...) where {D,T}
    fraction, initial_step = something(sys.compress)
    _check_compress(sys)
    o = _compress_center(sys)
    expansion = T(fraction)^(-one(T) / 3)
    _println()
    @_printf("  Packing at %.1f%% of the density (the geometry expanded by %.4f), then compressing\n",
        100 * fraction, expansion)
    # The constraints and the cell expanded (the molecules are placed by the packing)
    _scale_system!(sys, expansion, o; molecules=false)
    # The initial placement, at the low density
    sys.check_initial_placement && return _packmol(sys; run_options...)
    output_file = sys.output_file
    sys.output_file = ""
    scale = expansion  # the current scale of the geometry, relative to that of the input
    try
        _packmol(sys; run_options...)
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
        try
            while scale > one(T) + sqrt(eps(T))
                s = max(one(T) - step, one(T) / scale)
                saved_positions = copy(sys.molecule_positions)
                saved_constraints = [copy(st.constraints) for st in sys.structure_types]
                saved_unitcell = sys.unitcell
                _scale_system!(sys, s, o)
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
                else
                    # back to the last converged packing, with half the step
                    copyto!(sys.molecule_positions, saved_positions)
                    for (st, cs) in zip(sys.structure_types, saved_constraints)
                        copyto!(st.constraints, cs)
                    end
                    sys.unitcell = saved_unitcell
                    step /= 2
                    nconverged = 0
                    step < _compress_min_step && break
                end
            end
        finally
            _screen_log[] = screen_log
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
