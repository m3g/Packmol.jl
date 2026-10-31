#
# Move the worst molecules to new random positions.
# Scans fmol for free molecules with fmol > precision, computes the number
# of bad molecules and moves a fraction of them randomly. Returns the
# indices of the molecules actually moved, so the caller can act on them
# individually (e.g. packmol_main.jl fattens their atom radii back up).
# Following the Fortran Packmol heuristic (heuristics.f90 movebad subroutine).
#
# On top of that per-candidate probabilistic selection, the single
# worst-offending molecule of the system is additionally force-moved about
# half the time (one extra `rand()` draw, independent of the probabilistic
# scan below) — insurance against a molecule (or a pair of molecules tangled
# with each other) stuck at a persistent, but comparatively mild, violation
# that would otherwise keep drawing an arbitrarily small move probability
# (see the `fmol_max_type` comment below) call after call.
#
# Each moved molecule is placed by `_movebad_place_molecule!`: up to
# `max_guess_try` trials, each drawing a fresh random center of mass/rotation
# and then (when `fg_output`/`atom_positions`/`mol_iat_first` are supplied)
# running a short constraint-only fit on it, keeping whichever trial scores
# best (lowest post-fit constraint penalty among non-overlapping trials).
# This mirrors `reinitialize_with_bounds!`'s own best-of-N placement, with
# the constraint fit folded into each trial rather than left for a separate
# pass afterward.
#
function movebad!(
    packmol_system::PackmolSystem{D,T},
    fmol::Vector{T},
    free_mol_indices::Vector{Int},
    mol_structure_type::Vector{Int},
    RNG;
    movefrac::T=T(0.05),
    precision::T=T(1e-2),
    fmol_max_type::Vector{T}=zeros(T, length(packmol_system.structure_types)),
    cm_lo_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    cm_hi_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    fixed_sys = nothing,
    fixed_lo::SVector{D,T} = zero(SVector{D,T}),
    fixed_hi::SVector{D,T} = zero(SVector{D,T}),
    overlap_tol::T = packmol_system.tolerance,
    fg_output::Union{Nothing,InteratomicDistanceFG{D,T}} = nothing,
    atom_positions::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    mol_iat_first::Union{Nothing,Vector{Int}} = nothing,
    opt_nit::Int = 20,
    max_guess_try::Int = 20,
) where {D,T}
    ntypes = length(packmol_system.structure_types)
    bad_threshold = precision / packmol_system.nmols

    # The bad molecules of all types, and the worst one of each type and of
    # the whole system.
    candidates = Int[]
    fmol_max_now = zeros(T, ntypes)
    worst_imol = 0
    for imol in free_mol_indices
        fmol[imol] > bad_threshold || continue
        push!(candidates, imol)
        ist = mol_structure_type[imol]
        fmol_max_now[ist] = max(fmol_max_now[ist], fmol[imol])
        if worst_imol == 0 || fmol[imol] > fmol[worst_imol]
            worst_imol = imol
        end
    end
    nbad = length(candidates)
    nbad == 0 && return Int[]
    # `fmol_max_type` (owned and persisted by the caller across loops)
    # tracks the worst fmol ever observed for each type, not just this
    # call's own candidates: once the population of bad molecules of a type
    # becomes homogeneous late in the packing (all clustered near the same,
    # small fmol), using *this call's* max as the probability's reference
    # would make every candidate look nearly as bad as the worst one, driving
    # most of their move probabilities back up toward 0.5 even though none of
    # them is actually far from converged. Anchoring instead to the
    # historical worst-ever value keeps that ratio — and thus the
    # probability — small for a mildly-bad, tightly-clustered population,
    # exactly as it should be.
    fmol_max_type .= max.(fmol_max_type, fmol_max_now)

    # A single quota for the whole system: `movefrac` of the free molecules.
    # A quota per structure type (at least one molecule each) moves far more
    # molecules than that in systems of many small types (ten types of 6
    # molecules would move 10 molecules per call, instead of 3 if they were a
    # single type), disrupting the packing of the molecules that are almost
    # done, while a molecule stuck far from its region competes for its
    # type's quota only with its few siblings.
    nfree = length(free_mol_indices)
    nmove = max(1, min(nbad, round(Int, min(movefrac, nbad / nfree) * nfree)))

    # Each moved molecule is placed within the center-of-mass bounds of its type
    function place!(imol)
        ist = mol_structure_type[imol]
        lo, hi = if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
            l, h = cm_lo_type[ist], cm_hi_type[ist]
            all(l .< h) ? (l, h) : (nothing, nothing)
        else
            (nothing, nothing)
        end
        _movebad_place_molecule!(
            packmol_system, imol, packmol_system.structure_types[ist], RNG;
            cm_lo=lo, cm_hi=hi,
            fixed_sys, fixed_lo, fixed_hi, overlap_tol,
            fg_output, atom_positions, mol_structure_type, mol_iat_first,
            precision, opt_nit, max_guess_try,
        )
    end

    # On top of the probabilistic selection below — in which a persistently
    # stuck but comparatively mild offender can get an arbitrarily small move
    # probability, call after call — force-move the worst molecule of the
    # system about half the time. This is cheap insurance against exactly
    # that stuck case (e.g. two molecules genuinely tangled with each other,
    # unable to separate on their own) without making every call
    # deterministic.
    moved = Int[]
    forced_worst = rand(RNG, T) > T(0.5)
    if forced_worst
        place!(worst_imol)
        push!(moved, worst_imol)
    end

    # Move bad molecules randomly: probability of moving is proportional to
    # fmol, relative to the worst-ever fmol of the molecule's type (the worst
    # molecule of a type is moved with probability 0.5). The candidates of all
    # types are scanned in a shuffled order, not by molecule index: a scan
    # in index order that stops when the quota is filled would almost never
    # reach the types listed last when the first ones are large (e.g. 100
    # lipids after 1000 water molecules), leaving their molecules stuck
    # (confirmed on `bilayer-pbc.inp`: a lipid atom's constraint violation
    # stayed frozen over hundreds of loops). Shuffled, every type gets a share
    # of the moves proportional to its number of bad molecules.
    for imol in Random.shuffle(RNG, candidates)
        length(moved) >= nmove && break
        forced_worst && imol == worst_imol && continue
        prob = 0.5 * fmol[imol] / fmol_max_type[mol_structure_type[imol]]
        if rand(RNG, T) < prob
            place!(imol)
            push!(moved, imol)
        end
    end
    return moved
end

#
# Best-of-`max_guess_try` placement for a single molecule: each trial draws a
# random center of mass/rotation, optionally fits it against constraints only
# (when the scratch buffers are supplied), and scores it as the post-fit
# constraint penalty (or the raw one, when no fit is done) — with any trial
# that overlaps the fixed structure forced to the worst possible score,
# regardless of how good its own constraint penalty looks, since overlap
# with fixed atoms isn't reflected in that penalty at all. The best-scoring
# trial is kept; if every trial overlaps, the last trial is kept anyway
# (matching `randomize_molecule!`'s own fallback) rather than leaving the
# molecule at its original, already-bad position.
#
function _movebad_place_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    st::StructureType{D,T},
    RNG;
    cm_lo::Union{Nothing,SVector{D,T}},
    cm_hi::Union{Nothing,SVector{D,T}},
    fixed_sys,
    fixed_lo::SVector{D,T},
    fixed_hi::SVector{D,T},
    overlap_tol::T,
    fg_output::Union{Nothing,InteratomicDistanceFG{D,T}},
    atom_positions::Union{Nothing,Vector{SVector{D,T}}},
    mol_structure_type::Vector{Int},
    mol_iat_first::Union{Nothing,Vector{Int}},
    precision::T,
    opt_nit::Int,
    max_guess_try::Int,
) where {D,T}
    do_fit = !isnothing(fg_output) && !isnothing(atom_positions) && !isnothing(mol_iat_first)
    mol_list = [imol]
    overlap_positions = isnothing(fixed_sys) ? SVector{D,T}[] : Vector{SVector{D,T}}(undef, st.natoms)
    x = do_fit ? Vector{T}(undef, 2 * D) : T[]

    best_mp = packmol_system.molecule_positions[imol]
    best_score = typemax(T)
    last_mp = best_mp
    for _ in 1:max_guess_try
        mp = _random_molecule_position(packmol_system, RNG; cm_lo, cm_hi)
        packmol_system.molecule_positions[imol] = mp
        score = if do_fit
            x_mol = reinterpret(MoleculePosition{D,T}, x)
            x_mol[1] = mp
            spgresult = spgbox!(
                (g, x) -> constraint_only_fg_for_mols!(
                    g, x, fg_output, packmol_system, atom_positions, mol_list, mol_structure_type, mol_iat_first,
                ),
                x;
                nitmax=opt_nit,
                nfevalmax=10 * opt_nit,
                callback=(result) -> result.f < precision,
            )
            x_mol = reinterpret(MoleculePosition{D,T}, x)
            mp = x_mol[1]
            packmol_system.molecule_positions[imol] = mp
            spgresult.f
        else
            constraint_penalty_sum(mp, st)
        end
        last_mp = mp
        overlaps = !isnothing(fixed_sys) &&
            overlaps_fixed(mp, st.reference_coordinates, fixed_sys, overlap_positions, fixed_lo, fixed_hi, overlap_tol)
        score = overlaps ? typemax(T) : score
        if score < best_score
            best_score = score
            best_mp = mp
        end
        best_score < precision && break
    end
    packmol_system.molecule_positions[imol] = best_score < typemax(T) ? best_mp : last_mp
    return nothing
end
