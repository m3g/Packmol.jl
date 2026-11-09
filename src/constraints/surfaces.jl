#
# Surfaces given by points (with or without normals), as those of TS2CG or of a sampled
# mesh or curve, and the constraints of the distance to them:
#
#   within surface     <file> <distance> [offset]
#   not within surface <file> <distance> [offset]
#
# The distance is a signed distance `h` when the points have normals: an implicit moving
# least squares surface (Kolluri, ACM Trans. Algorithms 4, 18, 2008), the average of the
# distances to the tangent planes of the points within the smoothing radius `r`, weighted
# by (1 - d²/r²)³: smooth, positive on the side of the normals. Far from the surface (no
# point within `r`), it is the distance to the nearest point, signed by its normal, which
# is blended in with a small weight near the surface, for `h` to be continuous. Without
# normals, `h` is the distance to the nearest point (a curve, or a surface without sides).
#
# The points are in a grid of cells of side `r`, so that those within `r` of a position
# are in the 27 cells around it. Each cell also stores a point near it (a breadth-first
# propagation from the cells with points), from which the nearest point to a position far
# from the surface is found by a descent on the points.
#
export PointSurface, read_surface, surface_distance, WithinSurface, NotWithinSurface

weight_default[:surface] = 5.0

"""
    PointSurface(points, normals = nothing; radius = nothing)

A surface (or a curve) given by points, as an `AbstractVector` of 3D positions (Å),
optionally with the normals at the points (any length: they are normalized). With
normals, the distance to the surface is signed (positive on the side of the normals) and
smooth; without them, it is the distance to the nearest point.

`radius` is the smoothing radius of the signed distance: the points within this distance
of a position define the surface there. By default, 2.5 times the typical distance
between neighboring points.

See [`read_surface`](@ref), to read the points from a file, and [`WithinSurface`](@ref)
and [`NotWithinSurface`](@ref), for the constraints.
"""
struct PointSurface{T}
    points::Vector{SVector{3,T}}    # sorted by cell
    normals::Vector{SVector{3,T}}   # in the order of `points`; empty: no normals
    radius::T                       # smoothing radius, the side of the cells
    lo::SVector{3,T}                # lower corner of the grid
    ncells::NTuple{3,Int}
    cell_start::Vector{Int32}       # points of cell c: cell_start[c]:cell_start[c+1]-1
    seed::Vector{Int32}             # a point near each cell
    bbox::Tuple{SVector{3,T},SVector{3,T}}  # bounding box of the points
end

function PointSurface(points::AbstractVector, normals = nothing; radius = nothing)
    T = float(eltype(eltype(points)))
    return _point_surface(T, points, normals, radius)
end

function _point_surface(::Type{T}, points, normals, radius) where {T}
    pts = SVector{3,T}[SVector{3,T}(p[1], p[2], p[3]) for p in points]
    isempty(pts) && throw(ArgumentError("""\n
        A surface needs at least one point.

    """))
    nrm = SVector{3,T}[]
    if !isnothing(normals)
        length(normals) == length(pts) || throw(ArgumentError("""\n
            The number of normals ($(length(normals))) differs from the number of points ($(length(pts))).

        """))
        for n in normals
            v = SVector{3,T}(n[1], n[2], n[3])
            norm(v) > 0 || throw(ArgumentError("""\n
                A normal of the surface is zero.

            """))
            push!(nrm, v / norm(v))
        end
    end
    lo = reduce((a, b) -> min.(a, b), pts)
    hi = reduce((a, b) -> max.(a, b), pts)
    r = isnothing(radius) ? T(2.5) * _typical_spacing(pts, lo, hi) : T(radius)
    r > 0 || throw(ArgumentError("""\n
        The smoothing radius of a surface must be positive (got $r).

    """))
    glo = lo .- 2r
    ncells = Tuple(max.(1, ceil.(Int, (hi .+ 2r .- glo) ./ r)))
    prod(big.(ncells)) <= 500_000_000 || throw(ArgumentError("""\n
        The surface spans $(ncells) cells of side $r Å, too many: use a larger
        smoothing radius (`radius`).

    """))
    perm, cell_start = _cell_lists(pts, glo, r, ncells)
    pts = pts[perm]
    isempty(nrm) || (nrm = nrm[perm])
    seed = _seed_points(cell_start, ncells)
    return PointSurface{T}(pts, nrm, r, glo, ncells, cell_start, seed, (lo, hi))
end

@inline function _cell_of(lo::SVector{3}, side, ncells::NTuple{3,Int}, x)
    return ntuple(i -> clamp(floor(Int, (x[i] - lo[i]) / side) + 1, 1, ncells[i]), 3)
end
@inline _linear_cell(c::NTuple{3,Int}, n::NTuple{3,Int}) = c[1] + (c[2] - 1) * n[1] + (c[3] - 1) * n[1] * n[2]

# Points sorted by cell (the permutation), and the first point of each cell
function _cell_lists(pts::Vector{SVector{3,T}}, lo, side, ncells) where {T}
    cell = Int[_linear_cell(_cell_of(lo, side, ncells, p), ncells) for p in pts]
    perm = sortperm(cell)
    cell_start = zeros(Int32, prod(ncells) + 1)
    for c in cell
        cell_start[c+1] += 1
    end
    cell_start[1] = 1
    for c in 2:length(cell_start)
        cell_start[c] += cell_start[c-1]
    end
    return perm, cell_start
end

# Typical distance between neighboring points: the median distance to the nearest point,
# of a sample of the points, with cells of the size of the spacing of points filling the
# volume of their bounding box (an upper bound of that of a surface)
function _typical_spacing(pts::Vector{SVector{3,T}}, lo, hi) where {T}
    length(pts) < 2 && return one(T)
    ext = max.(hi .- lo, one(T))
    side = max(cbrt(prod(ext) / length(pts)), T(0.1))
    ncells = Tuple(max.(1, ceil.(Int, ext ./ side)))
    perm, cell_start = _cell_lists(pts, lo, side, ncells)
    sorted = pts[perm]
    step = max(1, length(pts) ÷ 2000)
    d = T[]
    for i in 1:step:length(sorted)
        x = sorted[i]
        best = typemax(T)
        for I in _near_cells(ncells, _cell_of(lo, side, ncells, x))
            lc = _linear_cell(Tuple(I), ncells)
            for m in cell_start[lc]:cell_start[lc+1]-1
                m != i && (best = min(best, sum(abs2, x - sorted[m])))
            end
        end
        best < typemax(T) && push!(d, sqrt(best))
    end
    isempty(d) && return side
    sort!(d)
    return d[(length(d)+1)÷2]
end

# The range of cells around the cell c (27, or less at the borders of the grid)
@inline function _near_cells(ncells::NTuple{3,Int}, c::NTuple{3,Int})
    return CartesianIndices((max(c[1] - 1, 1):min(c[1] + 1, ncells[1]),
        max(c[2] - 1, 1):min(c[2] + 1, ncells[2]), max(c[3] - 1, 1):min(c[3] + 1, ncells[3])))
end

# A point near each cell: breadth-first propagation, through the faces of the cells, of
# the first point of the cells that have points
function _seed_points(cell_start::Vector{Int32}, ncells::NTuple{3,Int})
    nc = prod(ncells)
    seed = zeros(Int32, nc)
    queue = Int[]
    for c in 1:nc
        if cell_start[c+1] > cell_start[c]
            seed[c] = cell_start[c]
            push!(queue, c)
        end
    end
    n1, n12 = ncells[1], ncells[1] * ncells[2]
    head = 1
    while head <= length(queue)
        c = queue[head]; head += 1
        i = (c - 1) % n1 + 1
        j = ((c - 1) ÷ n1) % ncells[2] + 1
        k = (c - 1) ÷ n12 + 1
        for (di, dj, dk) in ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))
            ii, jj, kk = i + di, j + dj, k + dk
            (1 <= ii <= ncells[1] && 1 <= jj <= ncells[2] && 1 <= kk <= ncells[3]) || continue
            cc = _linear_cell((ii, jj, kk), ncells)
            if seed[cc] == 0
                seed[cc] = seed[c]
                push!(queue, cc)
            end
        end
    end
    return seed
end

# The nearest point to x, by a descent from the point k (squared distance d2): moves to
# the nearest point to x among those of the cells around the current one, until none is
# nearer
function _descend(s::PointSurface, x, k::Int, d2)
    while true
        kk, dd = k, d2
        for I in _near_cells(s.ncells, _cell_of(s.lo, s.radius, s.ncells, s.points[k]))
            lc = _linear_cell(Tuple(I), s.ncells)
            for m in s.cell_start[lc]:s.cell_start[lc+1]-1
                e = sum(abs2, x - s.points[m])
                e < dd && (kk = m; dd = e)
            end
        end
        kk == k && return k, d2
        k, d2 = kk, dd
    end
end

#
# The distance `h` to the surface and its gradient. `near` is false when there is no
# point within the smoothing radius; then, if `skip_far` is true, the far search is
# skipped, and `h` is only known to be at least the radius (in absolute value).
#
function _surface_distance(s::PointSurface{T}, x::SVector{3}, skip_far::Bool = false) where {T}
    r2 = s.radius^2
    has_normals = !isempty(s.normals)
    A = zero(eltype(x)); B = zero(eltype(x))
    gA = zero(x); gB = zero(x)
    best = typemax(eltype(x)); k = 0
    for I in _near_cells(s.ncells, _cell_of(s.lo, s.radius, s.ncells, x))
        lc = _linear_cell(Tuple(I), s.ncells)
        for m in s.cell_start[lc]:s.cell_start[lc+1]-1
            v = x - s.points[m]
            d2 = sum(abs2, v)
            d2 < best && (best = d2; k = m)
            if has_normals && d2 < r2
                t = 1 - d2 / r2
                w = t^3
                gw = (-6 * t^2 / r2) * v
                n = s.normals[m]
                a = dot(n, v)
                A += w * a; B += w
                gA += a * gw + w * n
                gB += gw
            end
        end
    end
    near = best < r2
    if !near
        skip_far && return zero(eltype(x)), zero(x), false
        if k == 0
            k = Int(s.seed[_linear_cell(_cell_of(s.lo, s.radius, s.ncells, x), s.ncells)])
            best = sum(abs2, x - s.points[k])
        end
        k, best = _descend(s, x, k, best)
    end
    v = x - s.points[k]
    d = sqrt(best)
    if !has_normals
        return d, (d > 0 ? v / d : zero(x)), near
    end
    n = s.normals[k]
    sg = dot(n, v) >= 0 ? one(eltype(x)) : -one(eltype(x))
    hf = sg * d
    ghf = d > 0 ? sg * v / d : n
    # The distance to the nearest point, with a small weight: it takes over where the
    # weights of the points vanish (a crossfade over r/2 < d < r was no better in the
    # packings, and worse on a pinched mesh)
    ε = T(1e-3)
    A += ε * hf; B += ε; gA += ε * ghf
    h = A / B
    return h, (gA - h * gB) / B, near
end

"""
    surface_distance(s::PointSurface, x)

The distance of the position `x` to the surface `s`, signed if the surface has normals
(positive on their side), and its gradient.
"""
function surface_distance(s::PointSurface{T}, x) where {T}
    h, g, _ = _surface_distance(s, SVector{3,T}(x[1], x[2], x[3]))
    return h, g
end

function Base.show(io::IO, ::MIME"text/plain", s::PointSurface{T}) where {T}
    nrm = isempty(s.normals) ? "without normals" : "with normals"
    print(io, "PointSurface{$T}: $(length(s.points)) points $nrm, smoothing radius $(round(s.radius; digits = 2)) Å")
end
Base.show(io::IO, s::PointSurface{T}) where {T} = print(io, "PointSurface{$T}($(length(s.points)) points)")

"""
    read_surface(file; radius = nothing, T = Float64)

Reads the points of a surface (a [`PointSurface`](@ref)) from a file:

- a point file of [TS2CG](https://github.com/marrink-lab/TS2CG) (the `.dat` files written
  by its `PLM` module, as `point/OuterBM.dat`): positions and normals, in nm, converted to Å;
- a text file with a point per line: `x y z` or `x y z nx ny nz` (Å; lines starting with
  `#` are comments).

`radius` is the smoothing radius of the signed distance (see [`PointSurface`](@ref)).
"""
read_surface(file::AbstractString; radius = nothing, T::Type = Float64) = _read_surface(String(file), T, radius)

function _read_surface(file::String, ::Type{T}, radius) where {T}
    isfile(file) || throw(ArgumentError("""\n
        Surface file not found: $file

    """))
    pts = SVector{3,T}[]
    nrm = SVector{3,T}[]
    ts2cg = false
    nline = 0
    io = open(file)
    try
        while !eof(io)
            line = readline(io)
            nline += 1
            str = strip(line)
            (isempty(str) || startswith(str, "#")) && continue
            # TS2CG headers: "Box ...", "< Point NoPoints ...>", "< id domain_id area X Y Z Nx ... >"
            if startswith(str, "<") || startswith(str, "Box")
                ts2cg = true
                continue
            end
            words = split(str)
            vals = try
                [parse(T, w) for w in words]
            catch
                throw(ArgumentError("""\n
                    Could not read line $nline of the surface file $file:
                    $line

                """))
            end
            if ts2cg
                # id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2 (nm)
                length(vals) >= 9 || throw(ArgumentError("""\n
                    Line $nline of the TS2CG point file $file has $(length(vals)) values (expected 17).

                """))
                push!(pts, 10 * SVector{3,T}(vals[4], vals[5], vals[6]))
                push!(nrm, SVector{3,T}(vals[7], vals[8], vals[9]))
            elseif length(vals) == 3 || length(vals) == 6
                push!(pts, SVector{3,T}(vals[1], vals[2], vals[3]))
                if length(vals) == 6
                    push!(nrm, SVector{3,T}(vals[4], vals[5], vals[6]))
                end
            else
                throw(ArgumentError("""\n
                    Line $nline of the surface file $file has $(length(vals)) values:
                    expected 3 (x y z) or 6 (x y z nx ny nz).

                """))
            end
        end
    finally
        close(io)
    end
    if !isempty(nrm) && length(nrm) != length(pts)
        throw(ArgumentError("""\n
            Some points of the surface file $file have normals and others do not.

        """))
    end
    return _point_surface(T, pts, isempty(nrm) ? nothing : nrm, radius)
end

#
# The constraints
#
struct SurfaceDistance{Placement,T} <: Constraint{Placement,3,T}
    surface::PointSurface{T}
    distance::T
    offset::T
    weight::T
end

function SurfaceDistance{Placement}(surface::PointSurface{T}, distance, offset, weight) where {Placement,T}
    distance >= 0 || throw(ArgumentError("""\n
        The distance of a surface constraint must not be negative (got $distance).

    """))
    if isempty(surface.normals) && offset != 0
        throw(ArgumentError("""\n
            An offset of a surface constraint requires a surface with normals.

        """))
    end
    return SurfaceDistance{Placement,T}(surface, T(distance), T(offset), T(weight))
end

"""
    WithinSurface(surface::PointSurface, distance; offset = 0.0, weight = 5.0)

The atoms within `distance` (Å) of the surface, or of the surface parallel to it at the
signed distance `offset` along its normals: `|h - offset| ≤ distance`, where `h` is the
(signed) distance to the surface. Input file syntax:
`within surface file distance [offset]`.
"""
WithinSurface(surface::PointSurface, distance; offset = 0, weight = weight_default[:surface]) =
    SurfaceDistance{Inside}(surface, distance, offset, weight)

"""
    NotWithinSurface(surface::PointSurface, distance; offset = 0.0, weight = 5.0)

The atoms farther than `distance` (Å) from the surface, or from the surface parallel to
it at the signed distance `offset` along its normals: `|h - offset| ≥ distance`. Input
file syntax: `not within surface file distance [offset]`.
"""
NotWithinSurface(surface::PointSurface, distance; offset = 0, weight = weight_default[:surface]) =
    SurfaceDistance{Outside}(surface, distance, offset, weight)

function constraint_penalty(c::SurfaceDistance{Inside}, x)
    h, _, _ = _surface_distance(c.surface, x)
    v = abs(h - c.offset) - c.distance
    return v > 0 ? c.weight * v^2 : zero(v)
end

function constraint_gradient(c::SurfaceDistance{Inside}, x)
    h, g, _ = _surface_distance(c.surface, x)
    v = abs(h - c.offset) - c.distance
    return v > 0 ? (2 * c.weight * v * sign(h - c.offset)) * g : zero(x)
end

# Not within: no point within the smoothing radius means |h| ≥ radius, and then the
# constraint is satisfied if radius ≥ distance + |offset|: no need of the far search
_far_satisfies(c::SurfaceDistance) = c.surface.radius >= c.distance + abs(c.offset)

function constraint_penalty(c::SurfaceDistance{Outside}, x)
    h, _, near = _surface_distance(c.surface, x, _far_satisfies(c))
    (near || !_far_satisfies(c)) || return zero(eltype(x))
    v = c.distance - abs(h - c.offset)
    return v > 0 ? c.weight * v^2 : zero(v)
end

function constraint_gradient(c::SurfaceDistance{Outside}, x)
    h, g, near = _surface_distance(c.surface, x, _far_satisfies(c))
    (near || !_far_satisfies(c)) || return zero(x)
    v = c.distance - abs(h - c.offset)
    # sign(0) = 0: on the surface itself (offset 0) the direction is that of the normal
    s = h == c.offset ? one(v) : sign(h - c.offset)
    return v > 0 ? (-2 * c.weight * v * s) * g : zero(x)
end

function Base.show(io::IO, ::MIME"text/plain", c::SurfaceDistance{P,T}) where {P,T}
    name = P === Inside ? "WithinSurface" : "NotWithinSurface"
    what = P === Inside ? "within" : "not within"
    shift = iszero(c.offset) ? "" : " shifted by $(c.offset) Å along the normals"
    print(io, "$name{$T}: $what $(c.distance) Å of the surface of $(length(c.surface.points)) points$shift, weight = $(c.weight)")
end

#
# Input parsing: `within surface file distance [offset]` and `not within surface ...`,
# the file relative to the input file (the path is joined by the parser of the
# structure block). The surfaces read are shared by the constraints that name the same
# file.
#
const _surface_cache = Dict{Tuple{String,Float64},PointSurface{Float64}}()

function _cached_surface(file::String, ::Type{T}) where {T}
    T === Float64 || return _read_surface(file, T, nothing)
    isfile(file) || return _read_surface(file, T, nothing) # (throws)
    key = (abspath(file), mtime(file))
    return get!(() -> _read_surface(file, Float64, nothing), _surface_cache, key)
end

function _parse_surface_constraint(::Type{P}, keyword, data, ::Type{T}, filename) where {P,T}
    length(data) in (2, 3) || throw(ArgumentError("""\n
        Error parsing '$keyword' constraint data for $filename: expected
        `$keyword file distance [offset]`.

    """))
    distance, offset = try
        parse(T, data[2]), length(data) == 3 ? parse(T, data[3]) : zero(T)
    catch
        throw(ArgumentError("""\n
            Error parsing '$keyword' constraint data for $filename: expected
            `$keyword file distance [offset]`.

        """))
    end
    surface = _cached_surface(String(data[1]), T)
    return SurfaceDistance{P}(surface, distance, offset, T(weight_default[:surface]))
end

parse_constraint(::Val{Symbol("within surface")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T} =
    _parse_surface_constraint(Inside, "within surface", data, T, filename)
parse_constraint(::Val{Symbol("not within surface")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T} =
    _parse_surface_constraint(Outside, "not within surface", data, T, filename)

# The bounding box of the region of a `within` constraint (volume tools)
function _inside_box(c::SurfaceDistance{Inside})
    m = c.distance + abs(c.offset)
    return (c.surface.bbox[1] .- m, c.surface.bbox[2] .+ m)
end
