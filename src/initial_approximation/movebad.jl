#
# Move the worst molecules to new random positions.
# Scans fmol for free molecules with fmol > precision, computes the number
# of bad molecules and moves a fraction of them randomly. Returns the
# indices of the molecules actually moved, so the caller can act on them
# individually (e.g. packmol_main.jl fattens their atom radii back up).
# Following the Fortran Packmol heuristic (heuristics.f90 movebad subroutine).
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

    # Group free molecules by structure type once, so badness, quota and
    # selection are all computed *per type* below, rather than pooled across
    # the whole system. Pooling is what caused a real bug: with a single
    # system-wide quota (movefrac of the total free count) and one combined
    # scan that stops as soon as the quota is filled, a large structure type
    # (e.g. 1000 water molecules) fills that quota on its own almost every
    # call, so a much smaller type (e.g. 100 lipids) never gets scanned at
    # all — its worst molecules are never relocated and can get permanently
    # stuck at the same constraint violation. Confirmed by instrumentation on
    # `bilayer-pbc.inp`: one lipid atom's constraint violation stayed
    # bit-for-bit frozen across hundreds of outer loops. Giving each type its
    # own quota (movefrac of *that type's* free count) guarantees every type
    # gets a proportional share of moves every call, regardless of how the
    # other types are sized.
    free_by_type = [Int[] for _ in 1:ntypes]
    for imol in free_mol_indices
        push!(free_by_type[mol_structure_type[imol]], imol)
    end

    moved = Int[]
    for ist in 1:ntypes
        candidates = free_by_type[ist]
        nfree_type = length(candidates)
        nfree_type == 0 && continue

        nbad_type = 0
        fmol_max_now = zero(T)
        for imol in candidates
            if fmol[imol] > bad_threshold
                nbad_type += 1
                fmol_max_now = max(fmol_max_now, fmol[imol])
            end
        end
        nbad_type == 0 && continue
        # `fmol_max_type` (owned and persisted by the caller across loops)
        # tracks the worst fmol ever observed for this type, not just this
        # call's own candidates: once the population of bad molecules of this
        # type becomes homogeneous late in the packing (all clustered near
        # the same, small fmol), using *this call's* max as the probability's
        # reference would make every candidate look nearly as bad as the
        # worst one, driving most of their move probabilities back up toward
        # 0.5 even though none of them is actually far from converged.
        # Anchoring instead to the historical worst-ever value keeps that
        # ratio — and thus the probability — small for a mildly-bad,
        # tightly-clustered population, exactly as it should be.
        fmol_max_type[ist] = max(fmol_max_type[ist], fmol_max_now)

        # Number of molecules of this type to move
        frac = min(movefrac, nbad_type / nfree_type)
        nmove = max(1, min(nbad_type, round(Int, frac * nfree_type)))

        st = packmol_system.structure_types[ist]
        lo, hi = if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
            l, h = cm_lo_type[ist], cm_hi_type[ist]
            all(l .< h) ? (l, h) : (nothing, nothing)
        else
            (nothing, nothing)
        end

        # Move molecules of this type randomly: probability of moving is
        # proportional to fmol value (worse molecules are more likely to be
        # moved). Scanned in a shuffled order (not `candidates`' own
        # molecule-index order) so that, within this type, the early-exit
        # once `nmove` is filled doesn't systematically favor whichever
        # molecules happen to sit first.
        n_moved_type = 0
        for imol in Random.shuffle(RNG, candidates)
            n_moved_type >= nmove && break
            if fmol[imol] > bad_threshold
                # Probability increases with fmol value: move the worst-ever
                # molecule of this type with 0.5 probability, linearly
                # decreasing probability for better molecules.
                prob = 0.5 * fmol[imol] / fmol_max_type[ist]
                if rand(RNG, T) < prob
                    _movebad_place_molecule!(
                        packmol_system, imol, st, RNG;
                        cm_lo=lo, cm_hi=hi,
                        fixed_sys, fixed_lo, fixed_hi, overlap_tol,
                        fg_output, atom_positions, mol_structure_type, mol_iat_first,
                        precision, opt_nit, max_guess_try,
                    )
                    push!(moved, imol)
                    n_moved_type += 1
                end
            end
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

    # Hard bounds on the rotation-angle DOFs (constrain_rotation), matching
    # packmol_main.jl's own bounds on the main optimization: without these,
    # the constraint-only fit below is free to rotate the molecule away from
    # st.rotation_bounds while chasing a lower constraint penalty, silently
    # undoing the bounded draw from _random_molecule_position above.
    lower, upper = if do_fit && any(!isnothing, st.rotation_bounds)
        cm_lo_inf = SVector{D,T}(ntuple(_ -> T(-Inf), D))
        cm_hi_inf = SVector{D,T}(ntuple(_ -> T(Inf), D))
        ang_lo = SVector{D,T}(ntuple(d -> isnothing(st.rotation_bounds[d]) ? T(-Inf) : st.rotation_bounds[d][1], D))
        ang_hi = SVector{D,T}(ntuple(d -> isnothing(st.rotation_bounds[d]) ? T(Inf) : st.rotation_bounds[d][2], D))
        reinterpret(T, [MoleculePosition(cm_lo_inf, ang_lo)]), reinterpret(T, [MoleculePosition(cm_hi_inf, ang_hi)])
    else
        nothing, nothing
    end

    best_mp = packmol_system.molecule_positions[imol]
    best_score = typemax(T)
    last_mp = best_mp
    for _ in 1:max_guess_try
        mp = _random_molecule_position(packmol_system, RNG; cm_lo, cm_hi, rotation_bounds=st.rotation_bounds)
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
                lower, upper,
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
