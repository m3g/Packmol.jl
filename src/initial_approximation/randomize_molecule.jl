#
# Draw a single random candidate MoleculePosition (center of mass + Euler
# angles) within the placement region: the given per-type CM bounds when
# available (they're always tighter than — and, under PBC, often much
# smaller a fraction of — the full unit cell, since they come from a
# structure type's own constraints, e.g. a slab confined to a fraction of a
# much taller periodic box), otherwise the PBC unit cell when periodic,
# falling back to the large sidemax box when neither applies. Both callers
# already collapse degenerate bounds (extent ≤ 0 in some dimension) to
# `nothing` before calling this, so any `cm_lo`/`cm_hi` seen here is safe to
# use directly. The Euler angles are uniform in [0, 2π) (rotation constraints are
# penalties, satisfied by the constraint fit of `movebad!`). Shared by
# `randomize_molecule!` (single draw, accept/reject on fixed-atom overlap)
# and `movebad!`'s best-of-N placement (draw + constraint fit, keep the best
# scoring trial).
#
function _random_molecule_position(
    packmol_system::PackmolSystem{D,T},
    RNG;
    cm_lo::Union{Nothing,SVector{D,T}} = nothing,
    cm_hi::Union{Nothing,SVector{D,T}} = nothing,
) where {D,T}
    cm = if !isnothing(cm_lo) && !isnothing(cm_hi)
        has_pbc, uc, center = _typed_unitcell(packmol_system)
        _random_cm_in_cell(RNG, cm_lo, cm_hi - cm_lo, has_pbc, uc, center, packmol_system.wigner_seitz_cell)
    elseif !isnothing(packmol_system.unitcell)
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
        frac = SVector{D,T}(ntuple(_ -> rand(RNG, T) - T(0.5), D))
        SVector{D,T}(uc * frac) + center
    else
        sidemax = T(DEFAULT_SIDEMAX)
        SVector{D,T}(ntuple(_ -> sidemax * (T(2) * rand(RNG, T) - one(T)), D))
    end
    angles = _random_angles(RNG, SVector{D,T})
    return MoleculePosition(cm, angles)
end

#
# A uniform random center of mass in the box from `lo` to `lo + extent`. Under PBC, only
# the part of the box inside the cell centered at `center` (the Wigner-Seitz cell for
# dodecahedral and octahedral cells, see `_wrap_cm`) is sampled: a draw outside it is
# drawn again (up to `max_tries` draws). Taken as drawn, it would be wrapped into the
# cell later, adding its density to the region it maps to: the bounds of a type that
# fills the cell are its bounding box, and wrapping uniform draws of the bounding box of a
# truncated octahedron (46% of which fall outside it) gave 45% more water in some
# regions of the cell than in others.
#
function _random_cm_in_cell(
    rng, lo::SVector{D,T}, extent::SVector{D,T},
    has_pbc::Bool, uc::AbstractMatrix, center::SVector{D,T}, wigner_seitz::Bool;
    max_tries::Int = 100,
) where {D,T}
    cm = lo + SVector{D,T}(ntuple(d -> rand(rng, T) * extent[d], D))
    has_pbc || return cm
    for _ in 2:max_tries
        _wrap_cm(cm, uc, center, wigner_seitz) == cm && return cm
        cm = lo + SVector{D,T}(ntuple(d -> rand(rng, T) * extent[d], D))
    end
    return cm
end

#
# Randomly re-place a molecule within the placement region.
# Uses per-structure-type bounding box for non-PBC placement.
#
# If `fixed_sys` is given (see `_build_overlap_check_system`), up to
# `max_guess_try` random candidates are drawn and the first that doesn't
# come within `tol` of the fixed atoms is kept; if none do within budget,
# the last trial is kept anyway (matching the pre-existing single-draw
# behavior as a fallback, rather than looping indefinitely or leaving the
# molecule unmoved). `fixed_sys=nothing` (the default) skips the check
# entirely and draws once, exactly as before — `overlaps_fixed`'s
# `::Nothing` method always reports "no overlap", so no separate branch is
# needed here. `tol` should match whatever distance `fixed_sys` itself was
# built with (see `_build_overlap_check_system`'s caller).
#
function randomize_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    st::StructureType{D,T},
    RNG;
    cm_lo::Union{Nothing,SVector{D,T}} = nothing,
    cm_hi::Union{Nothing,SVector{D,T}} = nothing,
    fixed_sys = nothing,
    fixed_lo::SVector{D,T} = zero(SVector{D,T}),
    fixed_hi::SVector{D,T} = zero(SVector{D,T}),
    tol::T = packmol_system.tolerance,
    max_guess_try::Int = 20,
) where {D,T}
    return _randomize_molecule!(packmol_system, imol, st, RNG, cm_lo, cm_hi, fixed_sys, fixed_lo, fixed_hi, tol, max_guess_try)
end

# Positional version of `randomize_molecule!`: a keyword call whose arguments are not
# all of concrete types (e.g. `fixed_sys::Union{Nothing,ParticleSystem2}`) cannot be
# resolved at compile time, as `juliac --trim` requires.
function _randomize_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    st::StructureType{D,T},
    RNG,
    cm_lo::Union{Nothing,SVector{D,T}},
    cm_hi::Union{Nothing,SVector{D,T}},
    fixed_sys,
    fixed_lo::SVector{D,T},
    fixed_hi::SVector{D,T},
    tol::T,
    max_guess_try::Int=20,
) where {D,T}
    mol_positions = isnothing(fixed_sys) ? SVector{D,T}[] : Vector{SVector{D,T}}(undef, st.natoms)
    mp = MoleculePosition(zero(SVector{D,T}), zero(SVector{D,T}))
    for _ in 1:max_guess_try
        mp = _random_molecule_position(packmol_system, RNG; cm_lo, cm_hi)
        overlaps_fixed(mp, st.reference_coordinates, fixed_sys, mol_positions, fixed_lo, fixed_hi, tol) || break
    end
    packmol_system.molecule_positions[imol] = mp
    return nothing
end

#
# Random Euler angles of a uniformly distributed rotation. With the convention of
# `eulermat` (R = Rx(β) Ry(γ) Rz(θ)), β and θ are uniform in [0, 2π), and γ, in
# [-π/2, π/2], has density cos(γ)/2. Uniform angles (all three in [0, 2π)) are not
# uniform rotations: they leave the axes of the molecules depleted near the z axis.
#
_random_angles(rng, ::Type{SVector{3,T}}) where {T} =
    SVector{3,T}(T(2π) * rand(rng, T), asin(T(2) * rand(rng, T) - one(T)), T(2π) * rand(rng, T))
_random_angles(rng, ::Type{SVector{D,T}}) where {D,T} = SVector{D,T}(ntuple(_ -> T(2π) * rand(rng, T), D))
