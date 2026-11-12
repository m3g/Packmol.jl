#
# `voxel_place`: a direct port to Julia of the placement algorithm of Bentopy's `pack`
# step (Westendorp et al., Protein Science 2026, https://github.com/marrink-lab/bentopy,
# `src/pack/`), on Packmol's systems. The algorithm, and its effectiveness for crowded
# systems of large molecules, are Bentopy's; the port follows it as closely as Packmol's
# input allows:
#
# - The space is a grid of voxels of side `resolution` (the periodic cell, or the box of
#   the regions of the molecules), with a background of taken voxels: the fixed
#   molecules, and the molecules placed.
# - Each type has a compartment mask: the voxels whose corners satisfy the constraints
#   of the type that apply to all its atoms (Bentopy's masks; with the corners, the atoms
#   of a placed molecule satisfy the constraints, for convex regions). The voxels out of it are
#   taken for that type (Bentopy's session background).
# - A molecule is voxelized at a rotation conservatively, as `voxelize.rs`: the voxels
#   that intersect a sphere of the radius of each atom (Bentopy's bead radius; here the
#   radius of the atom, so that molecules that do not share voxels are at a distance
#   larger than the sum of the radii: the tolerance).
# - The types are placed in decreasing order of their moment (sum of the squared
#   distances of the atoms to their center, `--rearrange moment`), each type completely
#   before the next.
# - Each try puts the voxels of the molecule at a random free voxel of the session
#   background (the corner of the box of the molecule's voxels; without PBC, only where
#   the box fits in the space): the molecule is placed if none of its voxels is taken,
#   and its voxels are then taken. A new random rotation after each placement, and after
#   `max_tries ÷ max_tries_rot_div` unsuccessful tries with the same rotation. A type
#   gives up after `max_tries = max_tries_mult × number` unsuccessful tries, and the
#   molecules that were not placed are dropped.
#
# Adaptations to Packmol's input: constraints that apply to some atoms only, and
# rotation constraints, are checked on the atoms of a try that does not collide (a
# violation is an unsuccessful try; the rotations are drawn satisfying the rotation
# constraints); random free voxels are drawn uniformly (Bentopy pops them from a
# shuffled list of the free voxels, renewed when a tenth is used: the same
# distribution); only orthorhombic periodic cells.
#

struct VoxelSpace{T}
    origin::SVector{3,T}   # the lower corner
    h::SVector{3,T}        # the sides of a voxel
    n::NTuple{3,Int}
    periodic::Bool
end

# The voxels (indices from 0, relative to the lower corner of its box) of atoms at
# `x` with radii `r`: those that intersect a sphere of radius r around an atom.
# Returns the offsets and the corner of the box (in space coordinates, relative to
# the positions `x`).
function _voxelize(x::AbstractVector{SVector{3,T}}, r::AbstractVector{T}, h::SVector{3,T}) where {T}
    rmax = maximum(r)
    lo = reduce((a, b) -> min.(a, b), x) .- rmax
    hi = reduce((a, b) -> max.(a, b), x) .+ rmax
    dims = ntuple(k -> max(1, ceil(Int, (hi[k] - lo[k]) / h[k])), 3)
    vox = falses(dims)
    for (xj, rj) in zip(x, r)
        p = xj - lo
        l = ntuple(k -> max(0, floor(Int, (p[k] - rj) / h[k])), 3)
        u = ntuple(k -> min(dims[k] - 1, ceil(Int, (p[k] + rj) / h[k]) - 1), 3)
        r2 = rj^2
        for iz in l[3]:u[3], iy in l[2]:u[2], ix in l[1]:u[1]
            c = SVector{3,T}((ix + T(0.5)) * h[1], (iy + T(0.5)) * h[2], (iz + T(0.5)) * h[3])
            d = max.(abs.(p - c) .- h ./ 2, zero(T))
            sum(abs2, d) <= r2 && (vox[ix+1, iy+1, iz+1] = true)
        end
    end
    offsets = [(I[1] - 1, I[2] - 1, I[3] - 1) for I in findall(vox)]
    return offsets, dims, lo
end

# Linear index (1-based) of the voxel `p + o`, or 0 out of a non-periodic space
@inline function _voxel_linear(s::VoxelSpace, p::NTuple{3,Int}, o::NTuple{3,Int})
    i1, i2, i3 = p[1] + o[1], p[2] + o[2], p[3] + o[3]
    n1, n2, n3 = s.n
    if s.periodic
        i1, i2, i3 = mod(i1, n1), mod(i2, n2), mod(i3, n3)
    elseif !(0 <= i1 < n1 && 0 <= i2 < n2 && 0 <= i3 < n3)
        return 0
    end
    return 1 + i1 + n1 * (i2 + n2 * i3)
end

# Whether the voxels of a molecule (`offsets`, or `loffsets`, linear) at the voxel `p`
# fall on taken voxels. Without crossing the boundaries of the space (always so without
# PBC), the linear offsets: an addition per voxel.
@inline function _collides(bg::BitArray{3}, s::VoxelSpace, p::NTuple{3,Int}, dims, offsets, loffsets)
    n1, n2, n3 = s.n
    if p[1] + dims[1] <= n1 && p[2] + dims[2] <= n2 && p[3] + dims[3] <= n3
        base = 1 + p[1] + n1 * (p[2] + n2 * p[3])
        @inbounds for l in loffsets
            bg[base+l] && return true
        end
        return false
    end
    for o in offsets
        bg[_voxel_linear(s, p, o)] && return true
    end
    return false
end

# Take the voxels that intersect spheres of radii `r` at absolute positions `x` (the
# fixed molecules)
function _stamp_atoms!(bg::BitArray{3}, s::VoxelSpace{T}, x, r) where {T}
    for (xj, rj) in zip(x, r)
        p = xj - s.origin
        l = ntuple(k -> floor(Int, (p[k] - rj) / s.h[k]), 3)
        u = ntuple(k -> ceil(Int, (p[k] + rj) / s.h[k]) - 1, 3)
        for iz in l[3]:u[3], iy in l[2]:u[2], ix in l[1]:u[1]
            c = SVector{3,T}((ix + T(0.5)) * s.h[1], (iy + T(0.5)) * s.h[2], (iz + T(0.5)) * s.h[3])
            d = max.(abs.(p - c) .- s.h ./ 2, zero(T))
            sum(abs2, d) <= rj^2 || continue
            i = _voxel_linear(s, (ix, iy, iz), (0, 0, 0))
            i > 0 && (bg[i] = true)
        end
    end
    return bg
end

# The box of the region of a constraint (`nothing` if unbounded)
_constraint_extent(c) = nothing
_constraint_extent(c::Sphere{Inside}) = (c.center .- c.radius, c.center .+ c.radius)
_constraint_extent(c::Box{Inside}) = (c.center .- c.sides ./ 2, c.center .+ c.sides ./ 2)
_constraint_extent(c::Cube{Inside}) = (c.center .- c.side / 2, c.center .+ c.side / 2)
function _constraint_extent(c::Ellipsoid{Inside})
    e = abs(c.scale) .* SVector(c.a, c.b, c.c)
    return (c.center .- e, c.center .+ e)
end
function _constraint_extent(c::Cylinder{Inside})
    a = c.center + c.length * c.axis / norm(c.axis)
    return (min.(c.center, a) .- c.radius, max.(c.center, a) .+ c.radius)
end

# The constraints of a type that apply to all its atoms
function _whole_molecule_constraints(st)
    isempty(st.atom_constraints) && return collect(eachindex(st.constraints))
    return [ic for ic in eachindex(st.constraints) if all(l -> ic in l, st.atom_constraints)]
end

# The box of the region of a type (the intersection of the extents of its constraints)
function _type_extent(st::StructureType{D,T}, whole) where {D,T}
    lo = fill(T(-Inf), SVector{3,T}); hi = fill(T(Inf), SVector{3,T})
    for ic in whole
        e = _constraint_dispatch(_constraint_extent, st.constraints[ic])
        isnothing(e) && continue
        lo = max.(lo, e[1]); hi = min.(hi, e[2])
    end
    return lo, hi
end

# The compartment of the constraints `cs`: true for the voxels out of it
function _outside_mask(s::VoxelSpace{T}, cs, precision::T) where {T}
    out = trues(s.n)
    lo = fill(T(-Inf), SVector{3,T}); hi = fill(T(Inf), SVector{3,T})
    for c in cs
        e = _constraint_dispatch(_constraint_extent, c)
        isnothing(e) && continue
        lo = max.(lo, e[1]); hi = min.(hi, e[2])
    end
    rng = ntuple(3) do k
        a = isfinite(lo[k]) ? max(0, floor(Int, (lo[k] - s.origin[k]) / s.h[k])) : 0
        b = isfinite(hi[k]) ? min(s.n[k] - 1, ceil(Int, (hi[k] - s.origin[k]) / s.h[k])) : s.n[k] - 1
        a:b
    end
    # the corners of the voxels in the region (a voxel is in it if its 8 corners are)
    nc = map(length, rng) .+ 1
    corner_in = falses(nc)
    for kz in 0:nc[3]-1, ky in 0:nc[2]-1, kx in 0:nc[1]-1
        x = s.origin + SVector{3,T}((first(rng[1]) + kx) * s.h[1], (first(rng[2]) + ky) * s.h[2], (first(rng[3]) + kz) * s.h[3])
        f = zero(T)
        for c in cs
            f += _constraint_dispatch(constraint_penalty, c, x)
        end
        corner_in[kx+1, ky+1, kz+1] = f <= precision
    end
    for (kz, iz) in enumerate(rng[3]), (ky, iy) in enumerate(rng[2]), (kx, ix) in enumerate(rng[1])
        inside = true
        for dz in 0:1, dy in 0:1, dx in 0:1
            inside &= corner_in[kx+dx, ky+dy, kz+dz]
        end
        inside && (out[ix+1, iy+1, iz+1] = false)
    end
    return out
end

# A uniformly random rotation (a random unit quaternion), as Euler angles
function _random_rotation(RNG, ::Type{T}) where {T}
    q = SVector{4,T}(ntuple(_ -> randn(RNG, T), 4))
    w, x, y, z = q / norm(q)
    R = SMatrix{3,3,T,9}(
        1 - 2(y^2 + z^2), 2(x * y + w * z), 2(x * z - w * y),
        2(x * y - w * z), 1 - 2(x^2 + z^2), 2(y * z + w * x),
        2(x * z + w * y), 2(y * z - w * x), 1 - 2(x^2 + y^2),
    )
    return euler_angles(R)
end

#
# Place the free molecules of `packmol_system` (reference coordinates centered), Bentopy's
# way. Returns the positions of the placed molecules of each free type.
#
function _bentopy_pack(
    packmol_system::PackmolSystem{D,T}, RNG;
    resolution::T, radius::Union{Nothing,T}, max_tries_mult::Int, max_tries_rot_div::Int,
) where {D,T}
    sts = packmol_system.structure_types
    atoms = packmol_system.atoms
    free = [!st.fixed.fixed for st in sts]
    iat_first = cumsum([1; [st.number_of_molecules * st.natoms for st in sts]])[1:end-1]
    radii(ist) = isnothing(radius) ?
        [atoms[iat_first[ist]+j-1].radius for j in 1:sts[ist].natoms] : fill(radius, sts[ist].natoms)
    whole = [_whole_molecule_constraints(st) for st in sts]

    # The space: the periodic cell, or the box of the regions of the free types
    has_pbc, uc, center = _typed_unitcell(packmol_system)
    space = if has_pbc
        isdiag = all(uc[i, j] == 0 for i in 1:3, j in 1:3 if i != j)
        (isdiag && !packmol_system.wigner_seitz_cell) || throw(ArgumentError("""\n
            voxel_place: only orthorhombic periodic cells are supported.

        """))
        L = SVector{3,T}(uc[1, 1], uc[2, 2], uc[3, 3])
        n = ntuple(k -> max(1, round(Int, L[k] / resolution)), 3)
        VoxelSpace{T}(center - L / 2, L ./ SVector{3,T}(n), n, true)
    else
        lo = fill(T(Inf), SVector{3,T}); hi = fill(T(-Inf), SVector{3,T})
        for ist in eachindex(sts)
            free[ist] || continue
            l, u = _type_extent(sts[ist], whole[ist])
            (all(isfinite, l) && all(isfinite, u)) || throw(ArgumentError("""\n
                voxel_place: the region of the molecules of $(basename(sts[ist].filename)) is
                not bounded. Without periodic boundary conditions, each type needs a constraint
                that bounds it (inside box, cube, sphere, ellipsoid or cylinder), for the space
                of the grid.

            """))
            lo = min.(lo, l); hi = max.(hi, u)
        end
        n = ntuple(k -> max(1, ceil(Int, (hi[k] - lo[k]) / resolution)), 3)
        VoxelSpace{T}(lo, fill(resolution, SVector{3,T}), n, false)
    end
    @_printf("  Space: %d x %d x %d voxels of %.2f A\n", space.n..., resolution)

    # The global background: the fixed molecules
    background = falses(space.n)
    if packmol_system.avoid_overlap
        for ist in eachindex(sts)
            sts[ist].fixed.fixed || continue
            _stamp_atoms!(background, space, sts[ist].reference_coordinates, radii(ist))
        end
    end

    # The compartment masks, one per set of constraints
    precision = packmol_system.constraint_precision
    masks = Tuple{Vector{AnyConstraint{T}},BitArray{3}}[]
    mask_of = Dict{Int,Union{Nothing,BitArray{3}}}()
    for ist in eachindex(sts)
        free[ist] || continue
        cs = sts[ist].constraints[whole[ist]]
        if isempty(cs)
            mask_of[ist] = nothing
            continue
        end
        k = findfirst(m -> isequal(m[1], cs), masks)
        if isnothing(k)
            push!(masks, (cs, _outside_mask(space, cs, precision)))
            k = length(masks)
        end
        mask_of[ist] = masks[k][2]
    end
    @_printf("  Compartments: %d\n", length(masks))

    # The types by moment, the largest first (Bentopy's default `--rearrange moment`)
    moment(st) = sum(r -> sum(abs2, r), st.reference_coordinates)
    order_types = sort(findall(free); by = ist -> moment(sts[ist]), rev = true)
    session = similar(background)
    positions = Dict{Int,Vector{MoleculePosition{D,T}}}()
    nvox = prod(space.n)
    for ist in order_types
        st = sts[ist]
        ref = st.reference_coordinates
        r = radii(ist)
        m = mask_of[ist]
        isnothing(m) ? copyto!(session, background) : (session .= background .| m)
        partial = [ic for ic in eachindex(st.constraints) if !(ic in whole[ist])]
        target = st.number_of_molecules
        placed = MoleculePosition{D,T}[]
        max_tries = max_tries_mult * target
        max_tries_per_rotation = max(1, max_tries ÷ max_tries_rot_div)
        tries = 0
        tries_per_rotation = 0
        # a rotation, and the voxels of the molecule at it
        function new_rotation()
            angles = _random_rotation(RNG, T)
            for _ in 1:1000
                isempty(st.rotation_constraints) && break
                _rotation_penalty(st, eulermat(angles)) <= precision && break
                angles = _random_rotation(RNG, T)
            end
            R = eulermat(angles)
            x = [R * rj for rj in ref]
            offsets, dims, lo = _voxelize(x, r, space.h)
            loffsets = [o[1] + space.n[1] * (o[2] + space.n[2] * o[3]) for o in offsets]
            return angles, R, offsets, loffsets, dims, lo
        end
        angles, R, offsets, loffsets, dims, lo = new_rotation()
        while length(placed) < target
            tries >= max_tries && break
            if !space.periodic && any(k -> dims[k] > space.n[k], 1:3)
                tries += 1
                angles, R, offsets, loffsets, dims, lo = new_rotation()
                tries_per_rotation = 0
                continue
            end
            # a random free voxel (the ones taken are not tries)
            p = (0, 0, 0)
            found = false
            for _ in 1:10^7
                i = rand(RNG, 1:nvox) - 1
                session[i+1] && continue
                p = (i % space.n[1], (i ÷ space.n[1]) % space.n[2], i ÷ (space.n[1] * space.n[2]))
                if !space.periodic && (p[1] >= space.n[1] - dims[1] ||
                        p[2] >= space.n[2] - dims[2] || p[3] >= space.n[3] - dims[3])
                    continue
                end
                found = true
                break
            end
            found || break   # no free voxels left
            ok = !_collides(session, space, p, dims, offsets, loffsets)
            # the cm of the molecule at this position: the corner of its voxels at p
            cm = space.origin + SVector{3,T}(p) .* space.h - lo
            if ok && !isempty(partial)
                for (j, rj) in enumerate(ref)
                    isempty(st.atom_constraints) && break
                    x = R * rj + cm
                    f = zero(T)
                    for ic in st.atom_constraints[j]
                        ic in whole[ist] && continue
                        f += _constraint_dispatch(constraint_penalty, st.constraints[ic], x)
                    end
                    if f > precision
                        ok = false
                        break
                    end
                end
            end
            if !ok
                tries += 1
                tries_per_rotation += 1
                if tries_per_rotation >= max_tries_per_rotation
                    angles, R, offsets, loffsets, dims, lo = new_rotation()
                    tries_per_rotation = 0
                end
                continue
            end
            for o in offsets
                i = _voxel_linear(space, p, o)
                background[i] = true
                session[i] = true
            end
            push!(placed, MoleculePosition{D,T}(cm, angles))
            angles, R, offsets, loffsets, dims, lo = new_rotation()
            tries_per_rotation = 0
        end
        positions[ist] = placed
    end
    return positions
end

"""
    voxel_place(system::PackmolSystem; resolution, radius, max_tries_mult = 1000, max_tries_rot_div = 100, seed)
    voxel_place(input_file::String; kwargs...)

Place the free molecules of a system without overlaps, and drop the molecules that do
not fit: a direct port to Julia of the placement algorithm of
[Bentopy](https://github.com/marrink-lab/bentopy) (its packing step; Westendorp et al.,
*Protein Science* 2026, [doi:10.1002/pro.70480](https://doi.org/10.1002/pro.70480)),
whose algorithm this is.

The space is a grid of voxels of side `resolution` (default 2.5 times the largest atomic
radius of the free molecules; Bentopy's models use 5 Å voxels and beads of 2 Å).
The molecules are voxelized, the voxels that intersect a sphere of the radius of each
atom (or of `radius`, Bentopy's bead radius, if given), and placed by type, the types of
the largest moment first, each molecule at a random free voxel where none of its voxels
is taken by the fixed molecules or those already placed, and within the compartment of its
type: the voxels whose corners satisfy the constraints that apply to all its atoms. The
rotation changes after each placement and after `max_tries ÷ max_tries_rot_div`
unsuccessful tries; a type gives up after `max_tries = max_tries_mult × number`
unsuccessful tries. Without periodic boundary conditions, each type needs a constraint
that bounds its region (inside box, cube, sphere, ellipsoid or cylinder); periodic cells
must be orthorhombic.

Returns `(placed, dropped)`: `placed`, a new system with the molecules that fit, which
starts from their positions (its `restart_from` is set to them, so that
`packmol(placed)` starts from this placement), and `dropped`, the number of molecules of
each structure type of `system` that did not fit.

```julia
placed, dropped = voxel_place(system)
packmol(placed)
```
"""
function voxel_place(input_file::String; kwargs...)
    return voxel_place(read_packmol_input(input_file); kwargs...)
end

function voxel_place(
    system::PackmolSystem{D,T};
    resolution::Union{Nothing,Real} = nothing,
    radius::Union{Nothing,Real} = nothing,
    max_tries_mult::Integer = 1000,
    max_tries_rot_div::Integer = 100,
    seed::Integer = system.seed,
) where {D,T}
    D == 3 || throw(ArgumentError("""\n
        voxel_place: only 3D systems are supported.

    """))
    isnothing(system.restart_from) || throw(ArgumentError("""\n
        voxel_place: the system already starts from given positions (restart_from).

    """))
    t0 = time()
    RNG = seed < 0 ? Random.Xoshiro() : Random.Xoshiro(seed)
    sts = system.structure_types
    _center_reference_coordinates!(system)
    rmax = isnothing(radius) ?
        maximum(a.radius for a in system.atoms if !sts[a.structure_type_index].fixed.fixed) : T(radius)
    h = T(something(resolution, T(2.5) * rmax))
    _println("  Placing the molecules on a voxel grid (the algorithm of Bentopy)...")
    positions = _bentopy_pack(system, RNG;
        resolution = h, radius = isnothing(radius) ? nothing : T(radius),
        max_tries_mult = Int(max_tries_mult), max_tries_rot_div = Int(max_tries_rot_div))

    # The system of the placed molecules: the types with their numbers of placed
    # molecules (the types without any left out), starting from their positions
    dropped = zeros(Int, length(sts))
    kept_types = StructureType{D,T}[]
    kept_positions = MoleculePosition{D,T}[]
    imol = 0
    nfree = 0
    for (ist, st) in enumerate(sts)
        if st.fixed.fixed
            append!(kept_positions, system.molecule_positions[imol+1:imol+st.number_of_molecules])
            push!(kept_types, st)
        else
            nfree += st.number_of_molecules
            p = positions[ist]
            dropped[ist] = st.number_of_molecules - length(p)
            append!(kept_positions, p)
            isempty(p) || push!(kept_types, _with_number(st, length(p)))
        end
        imol += st.number_of_molecules
    end
    _renumber_conformations!(kept_types)
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(kept_types)
    copyto!(molecule_positions, kept_positions)
    placed = PackmolSystem{D,T}(;
        (f => getfield(system, f) for f in fieldnames(PackmolSystem{D,T})
            if !(f in (:structure_types, :atoms, :molecule_positions, :nmols, :restart_from, :status)))...,
        structure_types = kept_types, atoms, molecule_positions, nmols,
    )
    placed.restart_from = _restart_source(get_atoms(placed))

    ndropped = sum(dropped)
    @_printf("  Placed %d of %d free molecules in %.1f s", nfree - ndropped, nfree, time() - t0)
    _println(ndropped == 0 ? "." : "; did not fit:")
    for (ist, st) in enumerate(sts)
        dropped[ist] > 0 || continue
        @_printf("    %s: %d of %d\n", basename(st.filename), dropped[ist], st.number_of_molecules)
    end
    return placed, dropped
end

# A copy of a structure type with `n` molecules
function _with_number(st::StructureType{D,T}, n::Int) where {D,T}
    return StructureType{D,T}(;
        (f => getfield(st, f) for f in fieldnames(StructureType{D,T}) if f != :number_of_molecules)...,
        number_of_molecules = n,
    )
end

# The conformations of a molecule (consecutive types, `conformation = (i, n)`) renumbered
# after some were left out
function _renumber_conformations!(sts::Vector{StructureType{D,T}}) where {D,T}
    k = 1
    while k <= length(sts)
        if sts[k].conformation[2] == 1
            k += 1
            continue
        end
        # the set: consecutive types of increasing conformation index
        last = k
        while last < length(sts) && sts[last+1].conformation[2] > 1 &&
                sts[last+1].conformation[1] > sts[last].conformation[1]
            last += 1
        end
        n = last - k + 1
        for (i, j) in enumerate(k:last)
            sts[j] = StructureType{D,T}(;
                (f => getfield(sts[j], f) for f in fieldnames(StructureType{D,T}) if f != :conformation)...,
                conformation = (i, n),
            )
        end
        k = last + 1
    end
    return sts
end
