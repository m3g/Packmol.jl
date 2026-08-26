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
    # Count bad molecules and find fmol range among them
    nbad = 0
    fmol_max = zero(T)
    for imol in free_mol_indices
        if fmol[imol] > precision / packmol_system.nmols
            nbad += 1
            fmol_max = max(fmol_max, fmol[imol])
        end
    end
    nbad == 0 && return Int[]
    # Number of molecules to move
    frac = min(movefrac, nbad / nfree)
    nmove = max(1, min(nbad, round(Int, frac * nfree)))
    # Move molecules randomly: worse molecules are more likely to be moved
    # (see the probability formula below).
    moved = Int[]
    for imol in free_mol_indices
        length(moved) >= nmove && break
        if fmol[imol] > precision / packmol_system.nmols
            # Probability is `movefrac` (the target move fraction, e.g. 5%)
            # scaled down exponentially by how far this molecule's fmol falls
            # short of the worst one (fmol_max): in the degenerate case where
            # every bad molecule of this type is equally bad (fmol ==
            # fmol_max for all), the exponential factor is 1 for each of
            # them, so each independently has probability exactly `movefrac`
            # — the expected fraction moved is then exactly `movefrac`, never
            # more. In the normal case, where the worst molecule's fmol
            # genuinely stands out from the rest, the exponential pulls every
            # other candidate's probability well below `movefrac`, so the
            # expected fraction moved drops well under the target as the
            # population's badness becomes less uniform. The decay is
            # weighted by the molecule's own atom count (natoms_mol) so that
            # bigger molecules — whose fmol naturally accumulates more terms
            # simply from having more atoms — aren't penalized for that size
            # alone: the same absolute gap in fmol decays more gently for a
            # larger molecule.
            ist = mol_structure_type[imol]
            natoms_mol = packmol_system.structure_types[ist].natoms
            prob = movefrac * exp(-(fmol_max - fmol[imol]) / natoms_mol)
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
