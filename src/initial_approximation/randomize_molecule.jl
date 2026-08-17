#
# Randomly re-place a molecule within the placement region.
# Uses per-structure-type bounding box for non-PBC placement.
#
function randomize_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    st::StructureType{D,T},
    RNG;
    cm_lo::Union{Nothing,SVector{D,T}} = nothing,
    cm_hi::Union{Nothing,SVector{D,T}} = nothing,
) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
        frac = SVector{D,T}(ntuple(_ -> rand(RNG, T) - T(0.5), D))
        cm = SVector{D,T}(uc * frac) + center
    elseif !isnothing(cm_lo) && !isnothing(cm_hi)
        extent = cm_hi - cm_lo
        cm = cm_lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
    else
        sidemax = T(DEFAULT_SIDEMAX)
        cm = SVector{D,T}(ntuple(_ -> sidemax * (T(2) * rand(RNG, T) - one(T)), D))
    end
    angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(RNG, T), D))
    packmol_system.molecule_positions[imol] = MoleculePosition(cm, angles)
    return nothing
end
