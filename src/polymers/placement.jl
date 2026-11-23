#
# The placement of the polymers: initially, and when moved by movebad, each polymer is
# placed as a rigid molecule, in one of its conformations, and its monomers are put where
# they are in that conformation (aligned to its atoms). The rigid molecules are those of
# a system (`RigidPolymerSystem`) with the polymers as single molecules, and the other
# structure types as in the system.
#

struct RigidPolymerSystem{D,T}
    system::PackmolSystem{D,T}
    # The molecule of the system of each molecule of the rigid system (0 for a polymer)
    sys_molecule::Vector{Int}
    # The molecule of the first copy of each polymer, minus one
    polymer_offset::Vector{Int}
    # The structure type of each molecule of the rigid system
    mol_structure_type::Vector{Int}
end

function RigidPolymerSystem(sys::PackmolSystem{D,T}) where {D,T}
    polymer_of_type, _, _ = _type_monomers(sys)
    types = StructureType{D,T}[]
    sys_molecule = Int[]
    polymer_offset = zeros(Int, length(sys.polymers))
    imol = 0
    for (t, st) in enumerate(sys.structure_types)
        ip = polymer_of_type[t]
        if ip == 0
            push!(types, st)
            for i in 1:st.number_of_molecules
                push!(sys_molecule, imol + i)
            end
        elseif t == sys.polymers[ip].first_type
            polymer_offset[ip] = length(sys_molecule)
            for rst in sys.polymers[ip].rigid_types
                push!(types, rst)
                append!(sys_molecule, zeros(Int, rst.number_of_molecules))
            end
        end
        imol += st.number_of_molecules
    end
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(types)
    rsys = PackmolSystem{D,T}(;
        input_file=sys.input_file, output_file="", filetype=sys.filetype,
        tolerance=sys.tolerance, radscale=sys.radscale, structure_types=types,
        tolerance_precision=sys.tolerance_precision, constraint_precision=sys.constraint_precision,
        max_random_init=sys.max_random_init, adjust_constraints_on_init=sys.adjust_constraints_on_init,
        seed=sys.seed, avoid_overlap=sys.avoid_overlap, screen_log=sys.screen_log,
        unitcell=sys.unitcell, unitcell_center=sys.unitcell_center, wigner_seitz_cell=sys.wigner_seitz_cell,
        atoms, molecule_positions, nmols,
    )
    return RigidPolymerSystem{D,T}(rsys, sys_molecule, polymer_offset, _build_mol_structure_type(rsys))
end

#
# Puts the monomers of copy `q` of polymer `p` (index `ip`; copy `j` of its conformation
# `c`) where they are in its rigid molecule: each monomer aligned to its atoms
#
function _place_monomers!(
    sys::PackmolSystem{3,T}, rp::RigidPolymerSystem{3,T}, ip::Int, q::Int, c::Int, j::Int, mol_offset::Vector{Int},
) where {T}
    p = sys.polymers[ip]
    rmol = rp.polymer_offset[ip] + q
    rst = rp.system.structure_types[rp.mol_structure_type[rmol]]
    mp = rp.system.molecule_positions[rmol]
    R = eulermat(mp.angles)
    for m in 1:_nmonomers(p)
        st = sys.structure_types[_monomer_type(p, c, m)]
        observed = [R * rst.reference_coordinates[i] + mp.cm for i in p.monomer_atoms[m]]
        sys.molecule_positions[_monomer_molecule(p, mol_offset, c, m, j)] = _align_molecule(st.reference_coordinates, observed)
    end
    return sys
end

#
# The initial placement of a system with polymers: that of the rigid system, with the
# positions of its molecules copied to the system, and the monomers placed from the
# polymers. Returns whether the constraints are satisfied (see `set_initial_approximation!`).
#
function _polymer_initial_approximation!(sys::PackmolSystem{D,T}, rp::RigidPolymerSystem{D,T}, RNG; nloop::Int) where {D,T}
    rsys = rp.system
    _println("  Placing the polymers as rigid molecules")
    rfree = Int[]
    imol = 0
    for st in rsys.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            st.fixed.fixed || push!(rfree, imol)
        end
    end
    converged = set_initial_approximation!(rsys, rfree, RNG; nloop)
    # The rigid polymers packed before their monomers are released: polymers released
    # overlapping each other end up threaded through each other's monomers, and the
    # packing does not converge (a packing of the rigid polymers, converged or not, has
    # few overlaps). The packing need not converge: its warning is not shown
    if !sys.check_initial_placement
        _println("  Packing the polymers as rigid molecules")
        _packmol(rsys; restart=true, nloop, seed=sys.seed, warn_unconverged=false)
        rsys.status == :packing_ready ||
            _println("  The packing of the rigid polymers ended (not converged, as it may be)")
        _println("  Releasing the fragments of the flexible molecules")
    end
    _center_reference_coordinates!(sys)
    for (rmol, smol) in enumerate(rp.sys_molecule)
        smol > 0 && (sys.molecule_positions[smol] = rsys.molecule_positions[rmol])
    end
    mol_offset, _ = _type_offsets(sys)
    for (ip, p) in enumerate(sys.polymers)
        for (q, (c, j)) in enumerate(_polymer_copies(p))
            _place_monomers!(sys, rp, ip, q, c, j, mol_offset)
        end
    end
    return converged
end

#
# Movebad of the polymers: the polymers with monomers that do not satisfy the convergence
# criterion (`unconverged`) are shaken, each monomer displaced randomly by up to
# `displacement` and rotated by up to `rotation` (radians) about each axis, as the
# molecules of movebad: up to a fraction `movefrac` of the polymers, the worst one half
# of the times, each of the others with a probability proportional to its function.
# (Moving a polymer as a whole to a random position, in a dense system, put it across
# other polymers: the function increased tenfold.) Returns the number of polymers moved.
#
function _movebad_polymers!(
    sys::PackmolSystem{D,T}, fmol::Vector{T}, unconverged::Vector{Bool}, RNG;
    movefrac::T, displacement::T=sys.tolerance / 4, rotation::T=T(0.05),
) where {D,T}
    isempty(sys.polymers) && return 0
    mol_offset, _ = _type_offsets(sys)
    # The monomers of each bad polymer, and the function of the polymer
    bad = Vector{Int}[]
    fbad = T[]
    ntotal = 0
    for p in sys.polymers
        for (c, j) in _polymer_copies(p)
            ntotal += 1
            mols = [_monomer_molecule(p, mol_offset, c, m, j) for m in 1:_nmonomers(p)]
            if any(imol -> unconverged[imol], mols)
                push!(bad, mols)
                push!(fbad, sum(imol -> fmol[imol], mols))
            end
        end
    end
    isempty(bad) && return 0
    nmove = max(1, min(length(bad), round(Int, movefrac * ntotal)))
    fmax = maximum(fbad)
    function shake!(mols)
        for imol in mols
            mp = sys.molecule_positions[imol]
            dcm = displacement * (2 * rand(RNG, SVector{D,T}) .- 1)
            dangles = rotation * (2 * rand(RNG, SVector{D,T}) .- 1)
            sys.molecule_positions[imol] = MoleculePosition(mp.cm + dcm, mp.angles + dangles)
        end
        return nothing
    end
    moved = 0
    kworst = argmax(fbad)
    forced = rand(RNG, T) > T(0.5)
    if forced
        shake!(bad[kworst])
        moved += 1
    end
    for k in Random.shuffle(RNG, collect(eachindex(bad)))
        moved >= nmove && break
        forced && k == kworst && continue
        if rand(RNG, T) < T(0.5) * fbad[k] / fmax
            shake!(bad[k])
            moved += 1
        end
    end
    return moved
end

# The loops without a new minimum of the number of bad molecules after which the packing
# of polymers is stalled, and the number of trial placements of each stalled polymer
# (`_packmol`)
const _polymer_stall_loops = 20
const _polymer_stall_trials = 500

#
# The re-placement of stalled bad polymers (see `_packmol`): the worst bad polymers (a
# fraction `movefrac` of them, at least one) are placed again, whole, in the conformation
# of their structure file (no distortion of their past packing is carried), at the best
# of `ntrials` random positions fitted to their constraints: the one of least overlap with
# the atoms of all the other molecules, at their current positions. Their monomers are
# released from there.
#
# Why: the shake moves the monomers by a fraction of the tolerance, and a polymer trapped
# by its neighbors is not freed by it (in a vesicle of 5765 flexible POPC, a lipid whose
# chain end was held by its constraint, beyond its neighbors, kept a bond stretched by
# 0.66 Å for 300 loops). Tried: moving it whole to a random position fitted to its
# constraints, checking only the fixed molecules (that of movebad, in the system of the
# rigid polymers, whose positions are those of before the release), freed the lipid; its
# effect on the polymer-capped gold nanoparticle was unclear (1/3 runs, against 2/3 on
# main). Inflating the radii of the atoms of the stalled polymers (x2, decaying in 20
# loops), to push their neighbors away, did not free the lipid (16 stalls in 500 loops,
# the bond stretched by 0.68 Å throughout).
#
function _replace_stalled_polymers!(
    sys::PackmolSystem{D,T}, rp::RigidPolymerSystem{D,T}, fmol::Vector{T}, unconverged::Vector{Bool}, RNG;
    movefrac::T, precision::T, ntrials::Int=_polymer_stall_trials,
) where {D,T}
    isempty(sys.polymers) && return 0
    mol_offset, _ = _type_offsets(sys)
    # The bad polymers, the worst first: (function, polymer, copy, conformation, index)
    bad = Tuple{T,Int,Int,Int,Int}[]
    for (ip, p) in enumerate(sys.polymers)
        for (q, (c, j)) in enumerate(_polymer_copies(p))
            fp = zero(T)
            isbad = false
            for m in 1:_nmonomers(p)
                imol = _monomer_molecule(p, mol_offset, c, m, j)
                fp += fmol[imol]
                isbad |= unconverged[imol]
            end
            isbad && push!(bad, (fp, ip, q, c, j))
        end
    end
    isempty(bad) && return 0
    sort!(bad; rev=true)
    nmove = min(length(bad), max(1, round(Int, movefrac * length(bad))))
    # The fits to the constraints, in the system of the rigid polymers
    rsys = rp.system
    cm_lo, cm_hi = compute_cm_bounds(rsys)
    fg_output = InteratomicDistanceFG{D,T}(rsys)
    rsys_positions = Vector{SVector{D,T}}(undef, length(rsys.atoms))
    mol_iat_first = _build_mol_iat_first(rsys)
    # The atoms of the system, at their current positions
    atom_positions = Vector{SVector{D,T}}(undef, length(sys.atoms))
    topo = sys.polymer_topology
    for k in 1:nmove
        _, ip, q, c, j = bad[k]
        rmol = rp.polymer_offset[ip] + q
        ist = rp.mol_structure_type[rmol]
        st = rsys.structure_types[ist]
        instance = Int(topo.mol_instance[_monomer_molecule(sys.polymers[ip], mol_offset, c, 1, j)])
        current = _current_atoms(sys, atom_positions, rsys.structure_types)
        l, h = cm_lo[ist], cm_hi[ist]
        lo, hi = all(l .< h) ? _extended_cm_bounds(rsys, st, l, h) : (nothing, nothing)
        _movebad_place_molecule!(
            rsys, rmol, st, RNG;
            cm_lo=lo, cm_hi=hi, fixed_sys=nothing, fixed_lo=zero(SVector{D,T}), fixed_hi=zero(SVector{D,T}),
            overlap_tol=sys.tolerance, fg_output, atom_positions=rsys_positions, mol_structure_type=rp.mol_structure_type,
            mol_iat_first, precision, opt_nit=20, max_guess_try=ntrials, max_shift=T(Inf),
            current, skip_instance=instance, ntrials,
        )
        _place_monomers!(sys, rp, ip, q, c, j, mol_offset)
    end
    return nmove
end

# The options not supported with polymers
function _check_polymers(sys::PackmolSystem)
    isempty(sys.polymers) && return nothing
    isnothing(sys.restart_from) || throw(ArgumentError("""\n
        restart_from is not supported yet for systems with polymers.

    """))
    isnothing(sys.voxel_place) || throw(ArgumentError("""\n
        voxel_place is not supported for systems with polymers.

    """))
    return nothing
end
