#
# Guided movebad (experimental, off by default: `_guided_movebad[] = true`): the bad molecules
# are placed in the voids of the system instead of at random positions.
#
# A map of the voids (`VoidMap`): on a grid of spacing about `h` over the region of the atoms
# (the cell, for orthorhombic PBC), the distance `d(p)` of each point to the nearest atom,
# capped at the tolerance (only `max(0, tol - d)` matters), built by "splatting" each atom on
# the grid points within the tolerance. The molecules being moved are left out of it, and
# each one is added once placed (a local update), so that the next ones do not compete for
# the same void. The map is the same for all structure types: the cost of a trial position
# of a molecule is the overlap of its own atoms with the map, Σᵢ max(0, tol - d(xᵢ))², which
# accounts for its shape and orientation.
#
# Each trial draws the center of mass at a void point (`d ≥ tol`, within the bounds of the
# type) and the orientation at random, and is fitted to the constraints as in
# `_movebad_place_molecule!`; among the trials that satisfy the constraints, one is drawn
# with probability ∝ exp(-cost / T).
#
const _guided_movebad = Ref(false)
const _guided_spacing = Ref(1.0)       # grid spacing (Å)
const _guided_ntrials = Ref(20)        # trials per molecule
const _guided_temperature = Ref(0.1)   # T, relative to tol²

struct VoidMap{T}
    lo::SVector{3,T}
    h::SVector{3,T}
    n::NTuple{3,Int}
    periodic::Bool
    cap::T
    d::Array{T,3}
    voids::Vector{Int}     # linear indices of the grid points with d ≥ cap (at construction)
end

# The grid: the orthorhombic cell under PBC, else the bounding box of the atoms plus `cap`.
# `nothing` for other cells (triclinic, dodecahedral, octahedral: not implemented).
function _void_map(
    packmol_system::PackmolSystem{3,T}, atom_positions::Vector{SVector{3,T}},
    excluded_mols, mol_iat_first::Vector{Int}, h::T,
) where {T}
    cap = packmol_system.tolerance
    uc = packmol_system.unitcell
    periodic = !isnothing(uc)
    if periodic
        all(abs(uc[i, j]) < 1e-8 for i in 1:3, j in 1:3 if i != j) || return nothing
        L = SVector{3,T}(uc[1, 1], uc[2, 2], uc[3, 3])
        lo = packmol_system.unitcell_center - L / 2
        n = ntuple(k -> max(1, round(Int, L[k] / h)), 3)
        hs = SVector{3,T}(ntuple(k -> L[k] / n[k], 3))
    else
        lo = SVector{3,T}(typemax(T), typemax(T), typemax(T))
        hi = -lo
        for x in atom_positions
            lo = min.(lo, x); hi = max.(hi, x)
        end
        lo = lo .- cap; hi = hi .+ cap
        n = ntuple(k -> max(1, ceil(Int, (hi[k] - lo[k]) / h) + 1), 3)
        hs = SVector{3,T}(h, h, h)
    end
    d = fill(cap, n)
    vm = VoidMap{T}(lo, hs, n, periodic, cap, d, Int[])
    # atoms of the molecules not excluded
    excluded = falses(packmol_system.nmols)
    for imol in excluded_mols
        excluded[imol] = true
    end
    for iat in eachindex(atom_positions)
        excluded[packmol_system.atoms[iat].molecule_index] && continue
        _splat!(vm, atom_positions[iat])
    end
    for i in eachindex(d)
        d[i] >= cap && push!(vm.voids, i)
    end
    return vm
end

# The grid index of a coordinate (nearest point), wrapped under PBC; 0 out of the grid
@inline function _grid_index(vm::VoidMap, k::Int, x)
    i = round(Int, (x - vm.lo[k]) / vm.h[k]) + 1
    if vm.periodic
        return mod1(i, vm.n[k])
    else
        return 1 <= i <= vm.n[k] ? i : 0
    end
end
_grid_point(vm::VoidMap{T}, I::CartesianIndex{3}) where {T} =
    vm.lo + SVector{3,T}((I[1] - 1) * vm.h[1], (I[2] - 1) * vm.h[2], (I[3] - 1) * vm.h[3])

# Add an atom at `x`: the distances of the grid points within `cap` of it
function _splat!(vm::VoidMap{T}, x::SVector{3,T}) where {T}
    r = ntuple(k -> ceil(Int, vm.cap / vm.h[k]), 3)
    c = ntuple(k -> round(Int, (x[k] - vm.lo[k]) / vm.h[k]) + 1, 3)
    for i3 in c[3]-r[3]:c[3]+r[3], i2 in c[2]-r[2]:c[2]+r[2], i1 in c[1]-r[1]:c[1]+r[1]
        p = vm.lo + SVector{3,T}((i1 - 1) * vm.h[1], (i2 - 1) * vm.h[2], (i3 - 1) * vm.h[3])
        dist = norm(p - x)
        dist >= vm.cap && continue
        j1, j2, j3 = if vm.periodic
            mod1(i1, vm.n[1]), mod1(i2, vm.n[2]), mod1(i3, vm.n[3])
        else
            (1 <= i1 <= vm.n[1] && 1 <= i2 <= vm.n[2] && 1 <= i3 <= vm.n[3]) || continue
            i1, i2, i3
        end
        if dist < vm.d[j1, j2, j3]
            vm.d[j1, j2, j3] = dist
        end
    end
    return vm
end

# The distance to the nearest atom at `x` (the nearest grid point); `cap` out of the grid
function _void_distance(vm::VoidMap{T}, x::SVector{3,T}) where {T}
    i1 = _grid_index(vm, 1, x[1]); i2 = _grid_index(vm, 2, x[2]); i3 = _grid_index(vm, 3, x[3])
    (i1 == 0 || i2 == 0 || i3 == 0) && return vm.cap
    return vm.d[i1, i2, i3]
end

# The overlap cost of a molecule with the map
function _void_cost(vm::VoidMap{T}, mp::MoleculePosition{3,T}, st::StructureType{3,T}) where {T}
    R = eulermat(mp.angles)
    c = zero(T)
    for r in st.reference_coordinates
        δ = vm.cap - _void_distance(vm, R * r + mp.cm)
        δ > 0 && (c += δ^2)
    end
    return c
end

# A random void point within the bounds (up to `maxtries` draws), with a uniform jitter
# within its grid cell; `nothing` if none is found
function _random_void_cm(vm::VoidMap{T}, RNG, lo, hi; maxtries::Int = 1000) where {T}
    isempty(vm.voids) && return nothing
    CI = CartesianIndices(vm.d)
    for _ in 1:maxtries
        p = _grid_point(vm, CI[vm.voids[rand(RNG, 1:length(vm.voids))]])
        p = p + SVector{3,T}((rand(RNG, T) - T(0.5)) * vm.h[1], (rand(RNG, T) - T(0.5)) * vm.h[2], (rand(RNG, T) - T(0.5)) * vm.h[3])
        if isnothing(lo) || isnothing(hi) || all(lo .<= p .<= hi)
            return p
        end
    end
    return nothing
end

#
# Guided placement of a molecule (see the top of this file). Falls back to the random
# placement (`_movebad_place_molecule!`) when no void is found within the bounds.
#
function _guided_place_molecule!(
    vm::VoidMap{T},
    packmol_system::PackmolSystem{3,T},
    imol::Int,
    st::StructureType{3,T},
    RNG;
    cm_lo, cm_hi,
    fixed_sys, fixed_lo, fixed_hi, overlap_tol,
    fg_output, atom_positions, mol_structure_type, mol_iat_first,
    precision::T, opt_nit::Int, max_guess_try::Int, max_shift::T,
) where {T}
    mol_list = [imol]
    overlap_positions = isnothing(fixed_sys) ? SVector{3,T}[] : Vector{SVector{3,T}}(undef, st.natoms)
    x = Vector{T}(undef, 6)
    fit_lower, fit_upper = _fit_bounds(packmol_system, 1)
    trials = MoleculePosition{3,T}[]
    costs = T[]
    temperature = T(_guided_temperature[]) * packmol_system.tolerance^2
    for _ in 1:_guided_ntrials[]
        cm = _random_void_cm(vm, RNG, cm_lo, cm_hi)
        isnothing(cm) && break
        mp_drawn = MoleculePosition(cm, _random_angles(RNG, SVector{3,T}))
        x_mol = reinterpret(MoleculePosition{3,T}, x)
        x_mol[1] = mp_drawn
        packmol_system.molecule_positions[imol] = mp_drawn
        spgresult = spgbox!(
            (g, x) -> constraint_only_fg_for_mols!(
                g, x, fg_output, packmol_system, atom_positions, mol_list, mol_structure_type, mol_iat_first,
            ),
            x;
            lower=fit_lower, upper=fit_upper, nitmax=opt_nit, nfevalmax=10 * opt_nit,
            callback=(result) -> result.f < precision,
        )
        mp = reinterpret(MoleculePosition{3,T}, x)[1]
        spgresult.f < precision || continue
        norm(mp.cm - mp_drawn.cm) <= max_shift || continue
        !isnothing(fixed_sys) &&
            overlaps_fixed(mp, st.reference_coordinates, fixed_sys, overlap_positions, fixed_lo, fixed_hi, overlap_tol) &&
            continue
        push!(trials, mp)
        push!(costs, _void_cost(vm, mp, st))
    end
    if isempty(trials)
        _movebad_place_molecule!(
            packmol_system, imol, st, RNG;
            cm_lo, cm_hi, fixed_sys, fixed_lo, fixed_hi, overlap_tol,
            fg_output, atom_positions, mol_structure_type, mol_iat_first,
            precision, opt_nit, max_guess_try, max_shift,
        )
    else
        cmin = minimum(costs)
        w = [exp(-(c - cmin) / temperature) for c in costs]
        r = rand(RNG, T) * sum(w)
        k = 1
        while k < length(w) && r > w[k]
            r -= w[k]; k += 1
        end
        packmol_system.molecule_positions[imol] = trials[k]
    end
    # the molecule, in its new position, into the map
    mp = packmol_system.molecule_positions[imol]
    R = eulermat(mp.angles)
    for ref in st.reference_coordinates
        _splat!(vm, R * ref + mp.cm)
    end
    return nothing
end

#
# Guided initial placement (experimental, off by default: `_guided_init[] = true`), in place
# of the best-of-N random trials of `_place_molecules!`: the molecules are inserted one by
# one in a void map that starts with the fixed atoms, the structure types with more atoms
# first. Each trial draws the center of mass at a point of the map that is still empty
# (within the bounds of the type; uniform in them if none is found), and the orientation
# at random; trials that do not satisfy the constraints are fitted to them as in
# `_place_molecules!` (accepted if the fit moved the center of mass less than `max_shift`).
# The first trial that satisfies the constraints without overlapping the map is taken;
# else, after `_guided_ntrials[]` that satisfy them, one is drawn with probability
# ∝ exp(-cost / T). Each molecule placed is added to the map.
#
const _guided_init = Ref(false)

function _void_map_region(packmol_system::PackmolSystem{3,T}, lo::SVector{3,T}, hi::SVector{3,T}, h::T) where {T}
    cap = packmol_system.tolerance
    uc = packmol_system.unitcell
    if !isnothing(uc)
        all(abs(uc[i, j]) < 1e-8 for i in 1:3, j in 1:3 if i != j) || return nothing
        L = SVector{3,T}(uc[1, 1], uc[2, 2], uc[3, 3])
        lo = packmol_system.unitcell_center - L / 2
        n = ntuple(k -> max(1, round(Int, L[k] / h)), 3)
        hs = SVector{3,T}(ntuple(k -> L[k] / n[k], 3))
        return VoidMap{T}(lo, hs, n, true, cap, fill(cap, n), Int[])
    end
    lo = lo .- cap; hi = hi .+ cap
    n = ntuple(k -> max(1, ceil(Int, (hi[k] - lo[k]) / h) + 1), 3)
    return VoidMap{T}(lo, SVector{3,T}(h, h, h), n, false, cap, fill(cap, n), Int[])
end

function _guided_place_molecules!(
    packmol_system::PackmolSystem{3,T}, RNG, cm_min, cm_max, max_shift::T, opt_nit::Int, fit_precision::T,
    max_guess_try::Int, selected::Vector{Bool},
) where {T}
    sts = packmol_system.structure_types
    # bounds of each type, extended (as in `_place_molecules!`), and the region of the map
    has_pbc = !isnothing(packmol_system.unitcell)
    bounds = Vector{Tuple{SVector{3,T},SVector{3,T}}}(undef, length(sts))
    rlo = SVector{3,T}(typemax(T), typemax(T), typemax(T)); rhi = -rlo
    for (ist, st) in enumerate(sts)
        lo, hi = cm_min[ist], cm_max[ist]
        if all(hi - lo .> zero(T))
            lo, hi = _extended_cm_bounds(packmol_system, st, lo, hi)
        end
        bounds[ist] = (lo, hi)
        r = isempty(st.reference_coordinates) ? zero(T) : maximum(norm, st.reference_coordinates)
        rlo = min.(rlo, lo .- r); rhi = max.(rhi, hi .+ r)
    end
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{3,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    for iat in 1:natoms
        a = packmol_system.atoms[iat]
        if sts[a.structure_type_index].fixed.fixed
            rlo = min.(rlo, atom_positions[iat]); rhi = max.(rhi, atom_positions[iat])
        end
    end
    vm = _void_map_region(packmol_system, rlo, rhi, T(_guided_spacing[]))
    isnothing(vm) && return false
    for iat in 1:natoms
        a = packmol_system.atoms[iat]
        sts[a.structure_type_index].fixed.fixed && _splat!(vm, atom_positions[iat])
    end
    append!(vm.voids, 1:length(vm.d))

    _println("  Setting initial coordinates of the large molecules in the voids (guided, one at a time)...")
    progress_meter = ProgressBar(sum((st.number_of_molecules for (ist, st) in enumerate(sts) if selected[ist]); init = 0);
        desc="  Molecules placed: ")
    fg_output = InteratomicDistanceFG{3,T}(packmol_system)
    mol_structure_type = _build_mol_structure_type(packmol_system)
    mol_iat_first = _build_mol_iat_first(packmol_system)
    fit_bounds = _fit_bounds(packmol_system, 1)
    x = Vector{T}(undef, 6); mol_list = [0]
    temperature = T(_guided_temperature[]) * packmol_system.tolerance^2
    first_mol = cumsum([0; [st.number_of_molecules for st in sts]])
    order = sortperm([st.fixed.fixed ? 0 : st.natoms for st in sts]; rev = true)
    CI = CartesianIndices(vm.d)
    trials = MoleculePosition{3,T}[]; costs = T[]
    nzero = 0; nfeasible = 0
    for ist in order
        st = sts[ist]
        (st.fixed.fixed || !selected[ist]) && continue
        lo, hi = bounds[ist]
        for i in 1:st.number_of_molecules
            imol = first_mol[ist] + i
            empty!(trials); empty!(costs)
            best_mp = packmol_system.molecule_positions[imol]; best_fx = typemax(T)
            done = false
            for itry in 1:2 * max_guess_try
                # a center of mass in a void of the map, within the bounds
                cm = lo + SVector{3,T}(rand(RNG, T), rand(RNG, T), rand(RNG, T)) .* (hi - lo)
                for _ in 1:100
                    p = _grid_point(vm, CI[vm.voids[rand(RNG, 1:length(vm.voids))]])
                    p = p + SVector{3,T}((rand(RNG, T) - T(0.5)) * vm.h[1], (rand(RNG, T) - T(0.5)) * vm.h[2], (rand(RNG, T) - T(0.5)) * vm.h[3])
                    (all(lo .<= p .<= hi) && _void_distance(vm, p) >= vm.cap) || continue
                    cm = p; break
                end
                mp = MoleculePosition(cm, _random_angles(RNG, SVector{3,T}))
                fx = constraint_penalty_sum(mp, st, packmol_system)
                if fx >= fit_precision && itry > max_guess_try
                    fit_mp, fit_fx = _fit_to_constraints!(
                        x, mol_list, packmol_system, imol, mp, fg_output, atom_positions,
                        mol_structure_type, mol_iat_first, fit_precision, opt_nit, fit_bounds,
                    )
                    if fit_fx < fit_precision && norm(fit_mp.cm - mp.cm) <= max_shift
                        mp, fx = fit_mp, fit_fx
                    end
                end
                if fx < best_fx
                    best_fx = fx; best_mp = mp
                end
                fx < fit_precision || continue
                c = _void_cost(vm, mp, st)
                if c == 0
                    empty!(trials); empty!(costs); push!(trials, mp); push!(costs, c)
                    nzero += 1; done = true; break
                end
                push!(trials, mp); push!(costs, c)
                length(trials) >= _guided_ntrials[] && break
            end
            mp = if isempty(trials)
                best_mp
            elseif length(trials) == 1
                trials[1]
            else
                cmin = minimum(costs)
                w = [exp(-(c - cmin) / temperature) for c in costs]
                r = rand(RNG, T) * sum(w); k = 1
                while k < length(w) && r > w[k]
                    r -= w[k]; k += 1
                end
                trials[k]
            end
            isempty(trials) || (nfeasible += 1)
            packmol_system.molecule_positions[imol] = mp
            R = eulermat(mp.angles)
            for ref in st.reference_coordinates
                _splat!(vm, R * ref + mp.cm)
            end
            next!(progress_meter)
        end
    end
    finish!(progress_meter)
    _println("  Placed in voids: ", nzero, ", satisfying the constraints: ", nfeasible)
    return true
end

#
# Spreading of the centers of mass (experimental, `_spread_iters[] > 0`): after the trials of
# the initial placement, `_spread_iters[]` iterations of a short-ranged repulsion between the
# centers of mass of the free molecules, each a sphere of radius (tol/2) natoms^(1/3): a pair
# closer than the sum of their radii is pushed apart, each molecule by half the overlap
# (at most half its radius per iteration). The orientations are kept, and a move is
# accepted only if the molecule still satisfies its constraints (or does not violate them
# more than before). Orthorhombic cells and non-periodic systems only.
#
const _spread_iters = Ref(0)

function _spread_centers!(packmol_system::PackmolSystem{3,T}, free_mol_indices::Vector{Int}, precision::T) where {T}
    niter = _spread_iters[]
    niter > 0 || return nothing
    uc = packmol_system.unitcell
    L = if isnothing(uc)
        nothing
    else
        all(abs(uc[i, j]) < 1e-8 for i in 1:3, j in 1:3 if i != j) || return nothing
        SVector{3,T}(uc[1, 1], uc[2, 2], uc[3, 3])
    end
    sts = packmol_system.structure_types
    mol_structure_type = _build_mol_structure_type(packmol_system)
    radius = [packmol_system.tolerance / 2 * cbrt(T(st.natoms)) for st in sts]
    n = length(free_mol_indices)
    r = [radius[mol_structure_type[imol]] for imol in free_mol_indices]
    cutoff = 2 * maximum(r)
    disp = zeros(SVector{3,T}, n)
    t0 = time()
    accepted = 0; noverlap0 = 0; noverlap = 0
    for iter in 1:niter
        x = [packmol_system.molecule_positions[imol].cm for imol in free_mol_indices]
        pairs = isnothing(L) ? CellListMap.neighborlist(; positions = x, cutoff) :
            CellListMap.neighborlist(; positions = x, cutoff, unitcell = L)
        fill!(disp, zero(SVector{3,T}))
        nov = 0
        for (i, j, d) in pairs
            δ = r[i] + r[j] - d
            δ > 0 || continue
            nov += 1
            v = x[i] - x[j]
            isnothing(L) || (v = v - L .* round.(v ./ L))
            u = d > 0 ? v / d : SVector{3,T}(1, 0, 0)
            disp[i] += δ / 2 * u
            disp[j] -= δ / 2 * u
        end
        iter == 1 && (noverlap0 = nov)
        noverlap = nov
        nov == 0 && break
        acc = Threads.Atomic{Int}(0)
        Threads.@threads for k in 1:n
            dk = disp[k]
            nd = norm(dk)
            nd > 0 || continue
            nd > r[k] / 2 && (dk = dk * (r[k] / 2 / nd))
            imol = free_mol_indices[k]
            st = sts[mol_structure_type[imol]]
            mp = packmol_system.molecule_positions[imol]
            new = MoleculePosition(mp.cm + dk, mp.angles)
            p_old = constraint_penalty_sum(mp, st, packmol_system)
            p_new = constraint_penalty_sum(new, st, packmol_system)
            if p_new <= max(p_old, precision)
                packmol_system.molecule_positions[imol] = new
                Threads.atomic_add!(acc, 1)
            end
        end
        accepted += acc[]
    end
    @_printf("  Spread centers of mass: %d iterations, contacts %d -> %d, moves accepted %d, %.2f s\n",
        niter, noverlap0, noverlap, accepted, time() - t0)
    return nothing
end
