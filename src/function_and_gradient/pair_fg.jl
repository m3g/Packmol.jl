#
# The function and gradient of a pair of atoms closer than the sum of their radii,
# accumulated per molecule and per atom.
#
# Updates the function and gradient of the system given a pair of
# particles within the cutoff.
#
# `atom_radii` implements the same "loose start" heuristic as the original
# Fortran Packmol (getinp.f90's `discale`, default 1.1), generalized to a
# per-atom working radius rather than a single global scale factor: early in
# packing, every atom's own radius (thus its contribution to the pairwise
# target distance) is inflated, giving SPGBox a looser, easier target while
# molecules are far from feasible; each atom's entry in `atom_radii` decays
# toward its own floor (its user-specified or default radius) as improvement
# stalls, and is fattened back up when its molecule is relocated by
# movebad! (see packing/packing_loop.jl). `fg.dmin` tracks the true (unscaled)
# minimum distance and `fg.max_overlap` the largest overlap with the radii of the
# atoms, regardless of `atom_radii`, so convergence checks are unaffected by this
# internal optimization-target adjustment.
function cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system, atom_radii)
    (; x, y, i, j, d2) = pair
    iatom = packmol_system.atoms[i]
    jatom = packmol_system.atoms[j]
    if iatom.molecule_index == jatom.molecule_index
        return fg
    end
    tol = atom_radii[i] + atom_radii[j]
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    overlap = iatom.radius + jatom.radius - d
    fg.max_overlap = max(fg.max_overlap, overlap)
    if fg.count_overlaps[] && overlap >= fg.overlap_precision
        if overlap > fg.overlap_precision
            it, jt = minmax(iatom.structure_type_index, jatom.structure_type_index)
            fg.overlaps[it, jt] += 1
        end
        push!(fg.unconverged_mols, iatom.molecule_index, jatom.molecule_index)
    end
    if d < tol
        # Function value: add to total function value and to
        # molecular contributions
        fadd = (d - tol)^2
        fg.f += fadd
        _add_fmol!(fg, iatom.molecule_index, fadd)
        _add_fmol!(fg, jatom.molecule_index, fadd)
        # Gradient
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        _add_gxcar!(fg, i, -dvdd)
        _add_gxcar!(fg, j, dvdd)
    end
    return fg
end

# Accumulate into `fmol`/`gxcar`, recording the indices of the entries that
# become nonzero (see `touched_atoms`/`touched_mols`). An entry that returns to
# exactly zero and is touched again is listed twice, which is harmless: the
# reduction zeroes each entry of a copy after adding it to the main output.
@inline function _add_fmol!(fg::InteratomicDistanceFG, imol, v)
    fmol = fg.fmol[imol]
    iszero(fmol) && push!(fg.touched_mols, imol)
    fg.fmol[imol] = fmol + v
    return fg
end
@inline function _add_gxcar!(fg::InteratomicDistanceFG, iat, v)
    g = fg.gxcar[iat]
    iszero(g) && push!(fg.touched_atoms, iat)
    fg.gxcar[iat] = g + v
    return fg
end
