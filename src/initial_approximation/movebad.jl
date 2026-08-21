#
# Move the worst molecules to new random positions.
# Scans fmol for free molecules with fmol > precision, computes the number
# of bad molecules and moves a fraction of them randomly. Returns the
# indices of the molecules actually moved, so the caller can act on them
# individually (e.g. packmol_main.jl fattens their atom radii back up).
# Following the Fortran Packmol heuristic (heuristics.f90 movebad subroutine).
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
                    )
                else
                    randomize_molecule!(packmol_system, imol, st, RNG)
                end
                push!(moved, imol)
            end
        end
    end
    return moved
end