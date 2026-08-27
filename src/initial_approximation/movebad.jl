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
    loop::Int=0,
    nloop::Int=1,
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
    nfree = length(free_mol_indices)
    # Count bad molecules, the fmol range among them (per structure type,
    # both for this call alone and the historical worst-ever — see below),
    # and per-structure-type free-molecule totals (n_type) — the latter used
    # below to size each type's initial (loop 0) move quota.
    ntypes = length(packmol_system.structure_types)
    nbad = 0
    fmol_max_now = zeros(T, ntypes)
    n_type = zeros(Int, ntypes)
    for imol in free_mol_indices
        ist = mol_structure_type[imol]
        n_type[ist] += 1
        if fmol[imol] > precision / packmol_system.nmols
            nbad += 1
            fmol_max_now[ist] = max(fmol_max_now[ist], fmol[imol])
        end
    end
    nbad == 0 && return Int[]
    # `fmol_max_type` (owned and persisted by the caller across loops) tracks
    # the worst fmol ever observed for each structure type, not just this
    # call's own candidates: once the population of bad molecules becomes
    # homogeneous late in the packing (all clustered near the same, small
    # fmol), using *this call's* max as the exponential's reference would
    # make every candidate look nearly as bad as the worst one, driving most
    # of their probabilities back up toward the prefactor ceiling even
    # though none of them is actually far from converged. Anchoring instead
    # to the historical worst-ever value keeps the exponential gap — and
    # thus the probability — small for a mildly-bad, tightly-clustered
    # population, exactly as it should be.
    fmol_max_type .= max.(fmol_max_type, fmol_max_now)
    # Number of molecules to move
    frac = min(movefrac, nbad / nfree)
    nmove = max(1, min(nbad, round(Int, frac * nfree)))
    # As packing progresses, movebad! should disturb the system less: `stage`
    # ramps linearly from 0 (loop 0) to 1 (loop == nloop/2) and stays at 1
    # for the remainder of the run. It drives two things together, both
    # shrinking from their loop-0 behavior down to "at most 1 molecule of
    # each type" by the halfway point: the per-type move quota (quota_type,
    # a hard cap enforced below via moved_count_type) and the probability
    # prefactor itself (movefrac scaled by the same quota_type/quota0_type
    # ratio) — so a call late in the run is both less likely to move any
    # given candidate and structurally unable to move more than one member
    # of a type, rather than merely unlikely to.
    half = max(one(T), T(nloop) / 2)
    stage = clamp(T(loop) / half, zero(T), one(T))
    quota0_type = max.(1, round.(Int, movefrac .* n_type))
    quota_type = max.(1, round.(Int, (1 - stage) .* quota0_type .+ stage))
    moved_count_type = zeros(Int, ntypes)
    # Move molecules randomly: worse molecules are more likely to be moved
    # (see the probability formula below).
    moved = Int[]
    for imol in free_mol_indices
        length(moved) >= nmove && break
        if fmol[imol] > precision / packmol_system.nmols
            ist = mol_structure_type[imol]
            moved_count_type[ist] >= quota_type[ist] && continue
            # Probability is `movefrac` (the target move fraction, e.g. 5%),
            # itself shrunk by this call's stage-dependent quota ratio (see
            # above), scaled down exponentially by how far this molecule's
            # fmol falls short of the worst one ever seen for this type
            # (fmol_max_type[ist], the historical running max — see above,
            # not just this call's own candidates): in the degenerate case
            # where every bad molecule of this type is equally bad (fmol ==
            # fmol_max_type[ist] for all), the exponential factor is 1 for
            # each of them, so each independently has probability exactly
            # the (stage-scaled) prefactor — the expected fraction moved is
            # then exactly that prefactor, never more. In the normal case,
            # where the worst-ever fmol for this type genuinely stands out
            # from the current candidates, the exponential pulls every
            # candidate's probability well below the prefactor, so the
            # expected fraction moved drops well under the target as the
            # population's badness becomes less uniform (or simply smaller
            # than it once was). The decay is weighted by the molecule's own
            # atom count (natoms_mol) so that bigger molecules — whose fmol
            # naturally accumulates more terms simply from having more atoms
            # — aren't penalized for that size alone: the same absolute gap
            # in fmol decays more gently for a larger molecule.
            natoms_mol = packmol_system.structure_types[ist].natoms
            movefrac_now = movefrac * quota_type[ist] / quota0_type[ist]
            prob = movefrac_now * exp(-(fmol_max_type[ist] - fmol[imol]) / natoms_mol)
            if rand(RNG, T) < prob
                st = packmol_system.structure_types[ist]
                lo, hi = if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
                    l, h = cm_lo_type[ist], cm_hi_type[ist]
                    all(l .< h) ? (l, h) : (nothing, nothing)
                else
                    (nothing, nothing)
                end
                _movebad_place_molecule!(
                    packmol_system, imol, st, RNG;
                    cm_lo=lo, cm_hi=hi,
                    fixed_sys, fixed_lo, fixed_hi, overlap_tol,
                    fg_output, atom_positions, mol_structure_type, mol_iat_first,
                    precision, opt_nit, max_guess_try,
                )
                push!(moved, imol)
                moved_count_type[ist] += 1
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
