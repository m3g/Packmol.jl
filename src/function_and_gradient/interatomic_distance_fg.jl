using CellListMap: ParticleSystem, NeighborPair, pairwise!
using SPGBox: spgbox!

#
# Structure used to evaluate the function and gradient of the distance between atoms
# in a single pass, using CellListMap.
#
@kwdef mutable struct InteratomicDistanceFG{D,T}
    f::T
    g::Vector{MoleculePosition{D,T}}
    dmin::T
    fmol::Vector{T} # contribution of each molecule to the function value
    # Auxiliary array for gradient: carries the gradient relative to the cartesian coordinates of each atom
    gxcar::Vector{SVector{D,T}}
    # Pre-allocated per-thread buffers for fmol accumulation (avoids allocation in @spawn)
    fmol_threaded::Vector{Vector{T}}
    # Maximum constraint penalty over all atoms (for convergence checking)
    max_constraint_penalty::T = zero(T)
end

function InteratomicDistanceFG{D,T}(packmol_system::PackmolSystem) where {D,T}
    natoms = length(packmol_system.atoms)
    InteratomicDistanceFG(
        f = zero(T),
        g = zero(packmol_system.molecule_positions),
        dmin = typemax(T),
        fmol = zeros(T, packmol_system.nmols),
        gxcar = fill(zero(SVector{D,T}), natoms),
        fmol_threaded = [zeros(T, packmol_system.nmols) for _ in 1:Threads.nthreads()],
        max_constraint_penalty = zero(T),
    )
end

# Custom copy, reset and reducer functions for CellListMap.
# Note: `g` (molecule-level gradients) is excluded from copy/reduce — it is
# only written by chain_rule! after pairwise computation, so copying it per
# thread is wasted work (nmols × sizeof(MoleculePosition) per thread per call).
function CellListMap.copy_output(x::InteratomicDistanceFG{D,T}) where {D,T}
    InteratomicDistanceFG(
        x.f,
        MoleculePosition{D,T}[],  # empty — not used during pairwise
        x.dmin,
        copy(x.fmol),
        copy(x.gxcar),
        Vector{T}[],  # fmol_threaded — not used by CellListMap threads
        zero(T),  # max_constraint_penalty — not used by CellListMap threads
    )
end
function CellListMap.reset_output!(output::InteratomicDistanceFG{D,T}) where {D,T}
    output.f = zero(T)
    fill!(output.g, zero(MoleculePosition{D,T}))
    output.dmin = typemax(T)
    fill!(output.fmol, zero(T))
    fill!(output.gxcar, zero(SVector{D,T}))
    output.max_constraint_penalty = zero(T)
    return output
end
function CellListMap.reducer(x::InteratomicDistanceFG, y::InteratomicDistanceFG)
    x.f += y.f
    x.dmin = min(x.dmin, y.dmin)
    x.fmol .+= y.fmol
    x.gxcar .+= y.gxcar
    return x
end

# CellListMap's default reduce_output! (used by `reducer` above) merges the
# per-batch outputs with a plain serial `for ibatch in output_threaded`
# loop — for InteratomicDistanceFG that means a full O(nmols + natoms) vector
# add repeated once per batch, entirely on a single thread, on *every*
# pairwise! call (i.e. every fg! evaluation). At large scale (millions of
# atoms, tens of batches) this serial merge can itself take tens of
# milliseconds — measured at ~58ms for a 2.76M-atom/533K-molecule/10-batch
# system, i.e. a hard single-threaded tax paid on every optimizer iteration
# regardless of how many threads are available. This override instead
# parallelizes the merge by splitting the *index range* of fmol/gxcar across
# threads (each thread sums all batches' contributions for its own slice),
# which both uses all available threads and is more cache-friendly than the
# original batch-outer/index-inner loop order.
function CellListMap.reduce_output!(
    output::InteratomicDistanceFG{D,T}, output_threaded::Vector{InteratomicDistanceFG{D,T}},
) where {D,T}
    output.f = zero(T)
    output.dmin = typemax(T)
    for ot in output_threaded
        output.f += ot.f
        output.dmin = min(output.dmin, ot.dmin)
    end
    nmols = length(output.fmol)
    natoms = length(output.gxcar)
    @sync begin
        for (_, mrange) in enumerate(chunks(1:nmols; n=Threads.nthreads()))
            @spawn for i in mrange
                s = zero(T)
                for ot in output_threaded
                    s += ot.fmol[i]
                end
                output.fmol[i] = s
            end
        end
        for (_, arange) in enumerate(chunks(1:natoms; n=Threads.nthreads()))
            @spawn for i in arange
                s = zero(SVector{D,T})
                for ot in output_threaded
                    s += ot.gxcar[i]
                end
                output.gxcar[i] = s
            end
        end
    end
    return output
end

# Updates the function and gradient of the system given a pair of
# particles within the cutoff.
#
# `radscale` implements the same "loose start" heuristic as the original
# Fortran Packmol (getinp.f90's `discale`, default 1.1): early in packing,
# atom radii (thus the target distance) are inflated by `radscale`, giving
# GENCAN/SPGBox a looser, easier target while molecules are far from
# feasible; `radscale` decays toward 1.0 as improvement stalls (see
# packmol_main.jl). `fg.dmin` tracks the true (unscaled) minimum distance
# regardless of radscale, so convergence checks against the real tolerance
# are unaffected by this internal optimization-target adjustment.
function cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system, radscale)
    (; x, y, i, j, d2) = pair
    iatom = packmol_system.atoms[i]
    jatom = packmol_system.atoms[j]
    if iatom.molecule_index == jatom.molecule_index
        return fg
    end
    tol = radscale * (iatom.radius + jatom.radius)
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    if d < tol
        # Function value: add to total function value and to
        # molecular contributions
        fadd = (d - tol)^2
        fg.f += fadd
        fg.fmol[iatom.molecule_index] += fadd
        fg.fmol[jatom.molecule_index] += fadd
        # Gradient
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        fg.gxcar[i] -= dvdd
        fg.gxcar[j] += dvdd
    end
    return fg
end

#
# Compute Cartesian atomic positions from rigid-body molecule DOFs.
# atom_position[i] = R(angles) * reference_coord[i] + center_of_mass
#
function compute_atom_positions!(
    atom_positions::Vector{SVector{D,T}},
    molecule_positions,
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    imol_offset = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        ref = st.reference_coordinates
        natoms_st = st.natoms
        nmols_st = st.number_of_molecules
        st_imol_offset = imol_offset
        st_iat_offset = iat_offset
        @sync for (_, mol_range) in enumerate(chunks(1:nmols_st; n=Threads.nthreads()))
            @spawn begin
                for i in mol_range
                    imol = st_imol_offset + i
                    iat_first = st_iat_offset + (i - 1) * natoms_st
                    mp = molecule_positions[imol]
                    R = eulermat(mp.angles)
                    cm = mp.cm
                    for j in 1:natoms_st
                        atom_positions[iat_first + j] = R * ref[j] + cm
                    end
                end
            end
        end
        imol_offset += nmols_st
        iat_offset += nmols_st * natoms_st
    end
    return atom_positions
end

#
# Combined single-pass: compute Cartesian atom positions AND evaluate
# constraint penalties/gradients, for an explicit, possibly non-contiguous
# list of molecules (`mol_list`). Used by adjust_constraints!'s active-set
# loop: only molecules that are still in violation (or were just randomized)
# need their positions/constraints recomputed on a given pass — molecules
# already satisfying their constraints are untouched, since they are
# geometrically independent of everything else in this constraint-only phase.
#
# This does not accumulate into fg_output.f/.fmol/.gxcar via a whole-array
# reset + `+=`: since every atom in mol_list is written by exactly one task
# (mol_list is chunked into disjoint ranges) and constraint penalties/
# gradients depend only on that atom's own position, each entry is fully
# determined within this call — so it is assigned directly rather than
# accumulated, and no reset over molecules outside mol_list is needed. Only
# the shared scalars (fg_output.f, .max_constraint_penalty) need a lock
# across tasks.
function compute_positions_and_constraints_for_mols!(
    atom_positions::Vector{SVector{D,T}},
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
) where {D,T}
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    return _compute_positions_and_constraints_for_mols!(
        atom_positions, fg_output, packmol_system, mol_list, mol_structure_type, mol_iat_first,
        Val(has_pbc), unitcell, unitcell_center,
    )
end

# `HASPBC` is a type parameter (via Val), not a runtime Bool: `if HASPBC` below is
# resolved at compile time, so the non-PBC specialization never compiles the
# wrap_to_center branch at all. This matters because merely *having* that branch
# in the loop — even when never taken at runtime — makes the compiler infer `x`'s
# type as the join of both branches; since wrap_to_center's return type isn't
# perfectly inferred (it goes through a generic, allocating CellListMap path for
# a plain Matrix unitcell), that instability otherwise poisons every downstream
# constraint_penalty/constraint_gradient call in the loop (measured: ~650x more
# allocation for a mid-size, non-PBC system, on the full-system analog of this
# function — see _constraint_fg! below for the same effect measured directly).
function _compute_positions_and_constraints_for_mols!(
    atom_positions::Vector{SVector{D,T}},
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
) where {D,T,HASPBC}
    fg_output.f = zero(T)
    fg_output.max_constraint_penalty = zero(T)
    lk = ReentrantLock()
    @sync for (_, mrange) in enumerate(chunks(1:length(mol_list); n=Threads.nthreads()))
        @spawn begin
            f_local = zero(T)
            max_penalty_local = zero(T)
            for k in mrange
                imol = mol_list[k]
                ist = mol_structure_type[imol]
                st = packmol_system.structure_types[ist]
                ref = st.reference_coordinates
                natoms_st = st.natoms
                iat_first = mol_iat_first[imol]
                mp = packmol_system.molecule_positions[imol]
                R = eulermat(mp.angles)
                cm = mp.cm
                has_constraints = !isempty(st.constraints)
                fmol_local = zero(T)
                for j in 1:natoms_st
                    iat = iat_first + j - 1
                    pos = R * ref[j] + cm
                    atom_positions[iat] = pos
                    gxc = zero(SVector{D,T})
                    if has_constraints
                        x = HASPBC ? wrap_to_center(pos, unitcell, unitcell_center) : pos
                        atom = packmol_system.atoms[iat]
                        for ic in atom.constraints
                            c = st.constraints[ic]
                            penalty::T = constraint_penalty(c, x)::T
                            f_local += penalty
                            fmol_local += penalty
                            gxc += constraint_gradient(c, x)::SVector{D,T}
                        end
                    end
                    fg_output.gxcar[iat] = gxc
                end
                fg_output.fmol[imol] = fmol_local
                max_penalty_local = max(max_penalty_local, fmol_local)
            end
            @lock lk begin
                fg_output.f += f_local
                fg_output.max_constraint_penalty = max(fg_output.max_constraint_penalty, max_penalty_local)
            end
        end
    end
    return atom_positions
end

#
# Wrap position to the unit cell image centered at `center`.
# CellListMap.wrap_to_first maps to the cell with one corner at the origin;
# we shift so the cell is centered at `center` instead.
#
function wrap_to_center(x::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}) where {D,T}
    half = SVector{D,T}(ntuple(_ -> T(0.5), D))
    offset = SVector{D,T}(unitcell * half)
    return SVector{D,T}(CellListMap.wrap_to_first(x - center + offset, unitcell)) + center - offset
end

#
# Narrow `packmol_system.unitcell`/`unitcell_center` (typed `Union{Nothing,...}`,
# since PBC is optional) to concretely-typed local values *before* entering a
# @spawn loop that conditionally calls `wrap_to_center`. Reading the Union-typed
# field directly from inside a spawned closure — even guarded by an `if has_pbc`
# check computed outside it — prevents the compiler from proving the wrapped
# position stays a concrete SVector{D,T}, which cascades into heap-boxing every
# constraint_penalty/constraint_gradient call in the loop (measured: ~650x more
# allocation for a mid-size system). The dummy values are never read when
# `has_pbc` is false.
function _typed_unitcell(packmol_system::PackmolSystem{D,T}) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    unitcell::Matrix{T} = has_pbc ? packmol_system.unitcell : zeros(T, D, D)
    unitcell_center::SVector{D,T} = has_pbc ? packmol_system.unitcell_center : zero(SVector{D,T})
    return has_pbc, unitcell, unitcell_center
end

#
# Add constraint penalties and gradients to the fg output structure.
# When PBC is active, atom positions are wrapped to the unit cell centered
# at unitcell_center before evaluating constraints.
#
function constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    return _constraint_fg!(fg_output, atom_positions, packmol_system, Val(has_pbc), unitcell, unitcell_center)
end

# See the comment on _compute_positions_and_constraints_for_mols! above: HASPBC
# must be a type parameter (via Val), not a runtime Bool, so the non-PBC
# specialization never compiles the (allocating) wrap_to_center branch at all.
function _constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
) where {D,T,HASPBC}
    lk_fg = ReentrantLock()
    @sync for (ichunk, iat_range) in enumerate(index_chunks(packmol_system.atoms; n=Threads.nthreads()))
        @spawn begin
            fg_local = zero(fg_output.f)
            max_penalty_local = zero(T)
            fmol_local = fg_output.fmol_threaded[ichunk]
            fill!(fmol_local, zero(T))
            for iat in iat_range
                atom = packmol_system.atoms[iat]
                x = atom_positions[iat]
                if HASPBC
                    x = wrap_to_center(x, unitcell, unitcell_center)
                end
                st = packmol_system.structure_types[atom.structure_type_index]
                atom_penalty = zero(T)
                for ic in atom.constraints
                    c = st.constraints[ic]
                    penalty::T = constraint_penalty(c, x)::T
                    fg_local += penalty
                    fmol_local[atom.molecule_index] += penalty
                    atom_penalty += penalty
                    fg_output.gxcar[iat]::SVector{D,T} += constraint_gradient(c, x)::SVector{D,T}
                end
                max_penalty_local = max(max_penalty_local, atom_penalty)
            end
            @lock lk_fg begin
                fg_output.f += fg_local
                fg_output.fmol .+= fmol_local
                fg_output.max_constraint_penalty = max(fg_output.max_constraint_penalty, max_penalty_local)
            end
        end
    end
    return fg_output
end

#
# Combined objective function + gradient for SPGBox.
# Computes pairwise distance penalties + constraint penalties, then
# applies the chain rule to get gradients w.r.t. molecule DOFs.
#
function fg!(g, x,
    cl_system,
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    free_mol_indices::Vector{Int},
    radscale::T=one(T),
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    @sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=Threads.nthreads()))
        @spawn for k in krange
            packmol_system.molecule_positions[free_mol_indices[k]] = x_mol[k]
        end
    end
    # Compute Cartesian atomic coordinates from ALL molecule DOFs
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    # Update CellListMap positions
    @sync for (_, iatrange) in enumerate(chunks(eachindex(atom_positions); n=Threads.nthreads()))
        @spawn for iat in iatrange
            cl_system.xpositions[iat] = atom_positions[iat]
        end
    end
    # Compute pairwise distance penalties and Cartesian gradients
    pairwise!((pair, output) -> cartesian_fg!(pair, output, packmol_system, radscale), cl_system,)
    # Add constraint penalties and gradients
    constraint_fg!(cl_system.fg, atom_positions, packmol_system)
    # Zero molecule-level gradients (excluded from CellListMap reset/reduce)
    fill!(cl_system.fg.g, zero(MoleculePosition{D,T}))
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(cl_system.fg, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    @sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=Threads.nthreads()))
        @spawn for k in krange
            g_mol[k] = cl_system.fg.g[free_mol_indices[k]]
        end
    end
    return cl_system.fg.f
end

