#
# Move the worst molecules to new random positions.
# Scans fmol for free molecules with fmol > precision, computes the number
# of bad molecules and moves a fraction of them randomly. Returns the
# indices of the molecules actually moved, so the caller can act on them
# individually (e.g. packmol_main.jl fattens their atom radii back up).
# Following the Fortran Packmol heuristic (heuristics.f90 movebad subroutine).
#
# Each moved molecule goes through the same two-step placement as
# `adjust_constraints!`'s own movebad loop: (1) `randomize_molecule!` draws a
# new center of mass and rotation, rejecting (up to `max_guess_try` times)
# any trial that overlaps the fixed structure when `fixed_sys` is given; (2)
# once every molecule selected this call has been randomized, a short
# constraint-only optimization (`constraint_only_fg_for_mols!`, run only when
# `fg_output`/`atom_positions`/`mol_iat_first` are supplied) settles them into
# a constraint-satisfying position before control returns to the caller's own
# (distance-based) optimization. Molecules are geometrically independent in
# this constraint-only objective, so batching every moved molecule into one
# solver call is equivalent to — and cheaper than — minimizing them one at a
# time.
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
    # Move molecules randomly: probability of moving is proportional
    # to fmol value (worse molecules are more likely to be moved).
    moved = Int[]
    for imol in free_mol_indices
        length(moved) >= nmove && break
        if fmol[imol] > precision / packmol_system.nmols
            # Probability increases with fmol value: move the worst with 0.5
            # probablity, linearly decreasing probability for better molecules
            prob = 0.5 * fmol[imol] / fmol_max
            if rand(RNG, T) < prob
                ist = mol_structure_type[imol]
                st = packmol_system.structure_types[ist]
                if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
                    lo = cm_lo_type[ist]
                    hi = cm_hi_type[ist]
                    has_valid_bounds = all(lo .< hi)
                    randomize_molecule!(packmol_system, imol, st, RNG;
                        cm_lo = has_valid_bounds ? lo : nothing,
                        cm_hi = has_valid_bounds ? hi : nothing,
                        fixed_sys, fixed_lo, fixed_hi, tol = overlap_tol,
                    )
                else
                    randomize_molecule!(packmol_system, imol, st, RNG;
                        fixed_sys, fixed_lo, fixed_hi, tol = overlap_tol,
                    )
                end
                push!(moved, imol)
            end
        end
    end

    # Step 2: settle the just-randomized molecules into a constraint-proper
    # starting position (constraints only — this does not touch the
    # distance-based objective the caller optimizes afterward) before
    # returning control. Skipped when the caller didn't supply the scratch
    # buffers this needs (e.g. a system with no constraints at all wouldn't
    # benefit from it either).
    if !isempty(moved) && !isnothing(fg_output) && !isnothing(atom_positions) && !isnothing(mol_iat_first)
        x = Vector{T}(undef, length(moved) * 2 * D)
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(moved)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
        spgbox!(
            (g, x) -> constraint_only_fg_for_mols!(
                g, x, fg_output, packmol_system, atom_positions, moved, mol_structure_type, mol_iat_first,
            ),
            x;
            nitmax=opt_nit,
            nfevalmax=10 * opt_nit,
            callback=(result) -> result.f < precision,
        )
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(moved)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end
    end

    return moved
end