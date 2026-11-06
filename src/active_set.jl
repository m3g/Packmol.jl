#
# Packing loops on an active subset of the system.
#
# In a large system, after the first loops, most molecules are already packed: they
# contribute nothing to the objective function. The optimization of a loop can then be
# restricted to the molecules that are not (the "bad" molecules: overlaps or constraint
# violations) and to the molecules near them (any atom within `active_radius` of an atom
# of a bad molecule), with the other molecules frozen. Two cell lists: one of the active
# atoms, rebuilt at each evaluation, for the pairs of active atoms, and one of the frozen
# atoms (all the others, including those of the fixed molecules), built once per loop,
# crossed with the active atoms. Pairs of two frozen atoms are not computed (they do not
# change). Each evaluation of the objective function then costs in proportion to the
# active atoms, not to the whole system.
#
# The frozen molecules satisfy the tolerance and their constraints when the subsystem is
# built (otherwise they would be bad, thus active), and every pair of an active atom with
# a frozen one is computed: when the active molecules converge, the whole system does.
#

# Marks of molecules or atoms found within a distance of a set of atoms (cross cell list)
struct _Marks
    marks::Vector{Bool}
end
CellListMap.copy_output(m::_Marks) = _Marks(copy(m.marks))
CellListMap.reset_output!(m::_Marks) = (fill!(m.marks, false); m)
CellListMap.reducer(a::_Marks, b::_Marks) = (a.marks .|= b.marks; a)

# Marks, in `marks` (indexed by `index[j]` of the atom j of `ypositions`), every atom of
# `ypositions` within `cutoff` of an atom of `xpositions`
function _mark_within!(marks::Vector{Bool}, index::Vector{Int}, xpositions, ypositions, cutoff, unitcell)
    fill!(marks, false)
    isempty(xpositions) && return marks
    sys = ParticleSystem(
        xpositions=xpositions,
        ypositions=ypositions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=_Marks(marks),
        output_name=Val(:marks),
        parallel=_nbatches() > 1,
        nbatches=_cl_nbatches(),
        validate_coordinates=_validate_positions,
    )
    pairwise!((pair, m) -> (m.marks[index[pair.j]] = true; m), sys)
    marks .= sys.marks.marks
    return marks
end

mutable struct ActiveSubsystem{CS,CC}
    # Active (mobile) molecules, and their positions in the vector of free molecules
    mols::Vector{Int}
    kfree::Vector{Int}
    # Atoms of the active molecules (contiguous for each molecule, from `mol_first[k]`)
    atom_map::Vector{Int}
    mol_first::Vector{Int}
    nmobile::Int
    # Frozen atoms (the y set of the cross cell list)
    frozen_map::Vector{Int}
    cl::CS          # active atoms: pairs of active atoms, and the output of the subsystem
    cl_frozen::CC   # active atoms x frozen atoms (`nothing` if there are no frozen atoms)
end

function _active_output(packmol_system::PackmolSystem{D,T}, n::Int) where {D,T}
    return InteratomicDistanceFG{D,T}(
        f=zero(T),
        g=[zero(MoleculePosition{D,T})],  # not a thread copy (see `_is_thread_copy`)
        dmin=typemax(T),
        max_overlap=typemin(T),
        fmol=zeros(T, packmol_system.nmols),
        gxcar=fill(zero(SVector{D,T}), n),
        max_constraint_penalty=zero(T),
        overlap_precision=T(packmol_system.tolerance_precision),
    )
end

#
# Builds the active subsystem from the bad molecules (`unconverged`: those that do not
# satisfy the convergence criterion at the current positions) and the molecules near
# them, or `nothing` if the atoms of the active
# molecules are more than `max_fraction` of the atoms of the free molecules (the cost of
# an evaluation goes with the number of atoms; then the whole system is optimized), with
# the numbers of active molecules and of their atoms.
#
function _active_subsystem(
    packmol_system::PackmolSystem{D,T}, unconverged::Vector{Bool}, free_mol_indices::Vector{Int},
    cl_xpositions, unitcell, cutoff::T,
    mol_structure_type::Vector{Int}, mol_iat_first::Vector{Int};
    active_radius::T, max_fraction::T,
) where {D,T}
    atoms = packmol_system.atoms
    natoms = length(atoms)
    nmols = packmol_system.nmols
    cl_positions = SVector{D,T}[cl_xpositions[i] for i in 1:natoms]
    # The atoms of the bad molecules
    bad_atoms = SVector{D,T}[]
    nbad = 0
    for imol in free_mol_indices
        unconverged[imol] || continue
        nbad += 1
        st = packmol_system.structure_types[mol_structure_type[imol]]
        i0 = mol_iat_first[imol]
        for iat in i0:i0+st.natoms-1
            push!(bad_atoms, cl_positions[iat])
        end
    end
    nbad == 0 && return nothing, 0, 0
    # The molecules with any atom within active_radius of an atom of a bad molecule
    mol_of_atom = [a.molecule_index for a in atoms]
    active = Vector{Bool}(undef, nmols)
    _mark_within!(active, mol_of_atom, bad_atoms, cl_positions, active_radius, unitcell)
    is_free = falses(nmols)
    for imol in free_mol_indices
        is_free[imol] = true
    end
    mols = Int[]
    kfree = Int[]
    for (k, imol) in enumerate(free_mol_indices)
        if active[imol]
            push!(mols, imol)
            push!(kfree, k)
        end
    end
    nfree_atoms = 0
    nactive_atoms = 0
    for imol in free_mol_indices
        n = packmol_system.structure_types[mol_structure_type[imol]].natoms
        nfree_atoms += n
        active[imol] && (nactive_atoms += n)
    end
    nactive_atoms > max_fraction * nfree_atoms && return nothing, length(mols), nactive_atoms
    # Only free molecules move: the fixed ones are frozen even when near bad molecules
    for imol in 1:nmols
        is_free[imol] || (active[imol] = false)
    end
    atom_map = Int[]
    mol_first = Int[]
    for imol in mols
        st = packmol_system.structure_types[mol_structure_type[imol]]
        i0 = mol_iat_first[imol]
        push!(mol_first, length(atom_map) + 1)
        append!(atom_map, i0:i0+st.natoms-1)
    end
    nmobile = length(atom_map)
    frozen_map = [iat for iat in 1:natoms if !active[atoms[iat].molecule_index]]
    mobile_positions = [cl_positions[iat] for iat in atom_map]
    cl = ParticleSystem(
        xpositions=mobile_positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=_active_output(packmol_system, nmobile),
        output_name=Val(:fg),
        parallel=_nbatches() > 1,
        nbatches=_cl_nbatches(),
        validate_coordinates=_validate_positions,
    )
    cl_frozen = isempty(frozen_map) ? nothing : ParticleSystem(
        xpositions=copy(mobile_positions),
        ypositions=[cl_positions[iat] for iat in frozen_map],
        unitcell=unitcell,
        cutoff=cutoff,
        output=_active_output(packmol_system, nmobile),
        output_name=Val(:fg),
        parallel=_nbatches() > 1,
        nbatches=_cl_nbatches(),
        validate_coordinates=_validate_positions,
    )
    return ActiveSubsystem(mols, kfree, atom_map, mol_first, nmobile, frozen_map, cl, cl_frozen), length(mols), nmobile
end

# A pair of active atoms, as `cartesian_fg!`, with the local indices of the active atoms
function _active_cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system, atom_radii, sub::ActiveSubsystem)
    (; x, y, i, j, d2) = pair
    gi, gj = sub.atom_map[i], sub.atom_map[j]
    iatom = packmol_system.atoms[gi]
    jatom = packmol_system.atoms[gj]
    iatom.molecule_index == jatom.molecule_index && return fg
    tol = atom_radii[gi] + atom_radii[gj]
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    fg.max_overlap = max(fg.max_overlap, iatom.radius + jatom.radius - d)
    if d < tol
        fadd = (d - tol)^2
        fg.f += fadd
        _add_fmol!(fg, iatom.molecule_index, fadd)
        _add_fmol!(fg, jatom.molecule_index, fadd)
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        _add_gxcar!(fg, i, -dvdd)
        _add_gxcar!(fg, j, dvdd)
    end
    return fg
end

# A pair of an active atom (i, of x) with a frozen atom (j, of y): only the active atom
# has a gradient
function _frozen_cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system, atom_radii, sub::ActiveSubsystem)
    (; x, y, i, j, d2) = pair
    gi, gj = sub.atom_map[i], sub.frozen_map[j]
    iatom = packmol_system.atoms[gi]
    jatom = packmol_system.atoms[gj]
    tol = atom_radii[gi] + atom_radii[gj]
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    fg.max_overlap = max(fg.max_overlap, iatom.radius + jatom.radius - d)
    if d < tol
        fadd = (d - tol)^2
        fg.f += fadd
        _add_fmol!(fg, iatom.molecule_index, fadd)
        dv = y - x
        _add_gxcar!(fg, i, -2 * (d - tol) * dv / d)
    end
    return fg
end

#
# Objective function and gradient of the active subsystem, relative to the positions of
# the active molecules (`x`, in the order of `sub.mols`)
#
function _active_fg!(g, x, sub::ActiveSubsystem, packmol_system::PackmolSystem{D,T}, atom_radii::Vector{T}) where {D,T}
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    positions = sub.cl.xpositions
    # Positions of the atoms of the active molecules (as those of the cell list of the
    # whole system: each molecule shifted by the wrap of its center of mass)
    @_sync for krange in chunks(eachindex(sub.mols); n=_nbatches())
        @spawn for k in krange
            imol = sub.mols[k]
            mp = x_mol[k]
            packmol_system.molecule_positions[imol] = mp
            st = packmol_system.structure_types[packmol_system.atoms[sub.atom_map[sub.mol_first[k]]].structure_type_index]
            R = eulermat(mp.angles)
            shift = has_pbc ? _wrap_cm(mp.cm, unitcell, unitcell_center, packmol_system.wigner_seitz_cell) - mp.cm : zero(mp.cm)
            for j in 1:st.natoms
                positions[sub.mol_first[k] + j - 1] = R * st.reference_coordinates[j] + mp.cm + shift
            end
        end
    end
    pairwise!((pair, output) -> _active_cartesian_fg!(pair, output, packmol_system, atom_radii, sub), sub.cl)
    fg = sub.cl.fg
    # Pairs with the frozen atoms, added to the output of the subsystem
    if !isnothing(sub.cl_frozen)
        xf = sub.cl_frozen.xpositions
        for il in 1:sub.nmobile
            xf[il] = positions[il]
        end
        pairwise!((pair, output) -> _frozen_cartesian_fg!(pair, output, packmol_system, atom_radii, sub), sub.cl_frozen)
        ff = sub.cl_frozen.fg
        fg.f += ff.f
        fg.dmin = min(fg.dmin, ff.dmin)
        fg.max_overlap = max(fg.max_overlap, ff.max_overlap)
        for il in 1:sub.nmobile
            fg.gxcar[il] += ff.gxcar[il]
        end
    end
    # Constraints of the active atoms, and the chain rule
    lk = ReentrantLock()
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    @_sync for krange in chunks(eachindex(sub.mols); n=_nbatches())
        @spawn begin
            f_local = zero(T)
            max_local = zero(T)
            for k in krange
                imol = sub.mols[k]
                i0 = sub.mol_first[k]
                st = packmol_system.structure_types[packmol_system.atoms[sub.atom_map[i0]].structure_type_index]
                pmol = zero(T)
                for j in 1:st.natoms
                    il = i0 + j - 1
                    atom = packmol_system.atoms[sub.atom_map[il]]
                    xa = positions[il]
                    atom_penalty = zero(T)
                    for ic in atom.constraints
                        c = st.constraints[ic]
                        atom_penalty += _constraint_dispatch(constraint_penalty, c, xa)::T
                        fg.gxcar[il] += _constraint_dispatch(constraint_gradient, c, xa)::SVector{D,T}
                    end
                    pmol += atom_penalty
                    max_local = max(max_local, atom_penalty)
                end
                if !isempty(st.rotation_constraints)
                    p = _rotation_fg!(fg.gxcar, st, eulermat(x_mol[k].angles), i0)
                    pmol += p
                    max_local = max(max_local, p)
                end
                f_local += pmol
                g_mol[k] = partial_derivatives(
                    zero(MoleculePosition{D,T}), x_mol[k].angles, st.reference_coordinates,
                    @view(fg.gxcar[i0:i0+st.natoms-1]),
                )
            end
            @lock lk begin
                fg.f += f_local
                fg.max_constraint_penalty = max(fg.max_constraint_penalty, max_local)
            end
        end
    end
    return fg.f
end

# The rotation scaling of the active molecules (see rigid_body/rotation_metric.jl)
_subset(rs::RotationScaling{T}, k::Vector{Int}) where {T} =
    RotationScaling{T}(rs.inertia[k], rs.scaled[k], rs.Lt[k], rs.Lt_inv[k], rs.L_inv[k])
_subset(::Nothing, ::Vector{Int}) = nothing

