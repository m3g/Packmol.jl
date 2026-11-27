#
# Estimates for setting up a system: the volume that a set of molecules fills at a
# given density (`guess_volume`), the volume of the region available to the
# molecules of each structure type of an input, given its constraints and the fixed
# molecules (`free_volume`, by Monte Carlo), and the volume of the molecules themselves
# (`molecular_volume`, on a grid, as MoloVol).
#

const _AVOGADRO = 6.02214076e23

# Molar mass (g/mol) of a molecule
_molar_mass(atoms::AbstractVector{<:Atom}) = Float64(PDBToolsBase.mass(atoms))

"""
    guess_volume(input_file::String; density = 1.0)
    guess_volume(system::PackmolSystem; density = 1.0)
    guess_volume(molecules::AbstractVector{<:Pair}; density = 1.0)

The volume (Å³) that the molecules fill at the `density` (g/mL) of the whole system, and
the side (Å) of a cube of that volume, as the volume guesser of the Packmol web page. The
molecules are those of an input file or a `PackmolSystem` (fixed ones included), or
pairs of a structure file or a molar mass (g/mol) and a number of molecules:

```julia
guess_volume(["water.pdb" => 1000, "urea.pdb" => 400]; density = 1.1)
guess_volume([18.015 => 1000]) # 1000 water molecules at 1 g/mL
```

Returns `(; mass, volume, side)`: the total mass (g/mol, summed over the molecules), the
volume and the side of the cube, and prints a summary.
"""
function guess_volume(input_file::String; density::Real=1.0)
    return guess_volume(read_packmol_input(input_file, Val(3), Float64); density)
end

function guess_volume(system::PackmolSystem; density::Real=1.0)
    names = String[]
    masses = Float64[]
    numbers = Int[]
    for st in system.structure_types
        push!(names, basename(st.filename))
        push!(masses, _molar_mass(st.atoms))
        push!(numbers, st.number_of_molecules)
    end
    return _guess_volume(names, masses, numbers, Float64(density))
end

function guess_volume(molecules::AbstractVector{<:Pair}; density::Real=1.0)
    names = String[]
    masses = Float64[]
    numbers = Int[]
    for (molecule, number) in molecules
        if molecule isa AbstractString
            push!(names, basename(molecule))
            push!(masses, _molar_mass(read_pdb(expanduser(molecule))))
        else
            push!(names, string(molecule, " g/mol"))
            push!(masses, Float64(molecule))
        end
        push!(numbers, Int(number))
    end
    return _guess_volume(names, masses, numbers, Float64(density))
end

function _guess_volume(names::Vector{String}, masses::Vector{Float64}, numbers::Vector{Int}, density::Float64)
    density > 0 || throw(ArgumentError("\n\n    The density must be positive, got $density.\n\n"))
    mass = sum(masses[i] * numbers[i] for i in eachindex(masses); init=0.0)
    # g/mol -> g, then mL -> Å³ (1 mL = 1e24 Å³)
    volume = mass / _AVOGADRO / density * 1e24
    side = cbrt(volume)
    _println()
    _println("  Volume guess, at a density of ", _fmt(density, 4), " g/mL:")
    _println()
    for i in eachindex(names)
        _println("    ", _rpad(names[i], 30), _lpad(string(numbers[i]), 9), " molecules of ", _fmt(masses[i], 2), " g/mol")
    end
    _println()
    _println("  Total mass:   ", _fmt(mass, 2), " g/mol")
    _println("  Volume:       ", _fmt(volume, 1), " Å³")
    _println("  Side of cube: ", _fmt(side, 2), " Å")
    _println()
    return (; mass, volume, side)
end

# `s` padded with spaces to `n` characters, on the left or on the right (`lpad` and
# `rpad` don't compile with `juliac --trim`)
function _pad(s::AbstractString, n::Int, left::Bool)
    out = String(s)
    for _ in 1:(n - length(out))
        out = left ? " " * out : out * " "
    end
    return out
end
_lpad(s::AbstractString, n::Int) = _pad(s, n, true)
_rpad(s::AbstractString, n::Int) = _pad(s, n, false)

# A number with `digits` decimals, for the summaries (an integer from 10⁴ on)
_fmt(x::Real, digits::Int) = !isfinite(x) ? string(x) : abs(x) >= 1e4 ? string(round(Int, x)) : string(round(x; digits))

"""
    free_volume(input_file::String; npoints = 1_000_000, seed = 1)
    free_volume(system::PackmolSystem; npoints = 1_000_000, seed = 1)

The volume of the region available to each structure type, by Monte Carlo.

The volume (Å³) of the region of the molecules of each structure type of an input file or
a `PackmolSystem` is estimated from `npoints` random points of a box that contains the
region: the fraction of them that satisfies the constraints of the structure type and is
farther than the tolerance from the atoms of the fixed molecules.
With the volume, the density at which the molecules of the structure type would fill the
region (g/mL), which shows when too many (or too few) molecules were asked for.

The region is that of the constraints that apply to whole molecules. Those of some atoms
only (in `atoms` blocks) are ignored: the volume is then an upper bound, and a structure
type with constraints of some atoms only (as the lipids of a membrane, with their heads
and tails above and below planes) has no estimate. Neither has one a structure type whose
region is not bounded, by the periodic cell or by an `inside` box, cube, sphere,
cylinder or ellipsoid.

Returns a vector with, for each structure type, `(; structure, number, volume, error,
density, note)`: the error is the standard error of the estimate (Å³), and `volume`,
`error` and `density` are `NaN` without an estimate (`note` says why). Prints a summary.
"""
function free_volume(input_file::String; npoints::Integer=1_000_000, seed::Integer=1)
    return free_volume(read_packmol_input(input_file, Val(3), Float64); npoints, seed)
end

function free_volume(system::PackmolSystem{3,T}; npoints::Integer=1_000_000, seed::Integer=1) where {T}
    npoints > 0 || throw(ArgumentError("\n\n    npoints must be positive, got $npoints.\n\n"))
    rng = Random.Xoshiro(seed)
    has_pbc, uc, center = _typed_unitcell(system)
    fixed = SVector{3,T}[]
    for st in system.structure_types
        st.fixed.fixed && append!(fixed, st.reference_coordinates)
    end
    tol = system.tolerance
    results = NamedTuple{(:structure, :number, :volume, :error, :density, :note),Tuple{String,Int,Float64,Float64,Float64,String}}[]
    for st in system.structure_types
        st.fixed.fixed && continue
        name = basename(st.filename)
        # The constraints of the whole molecule: those of every atom
        whole = isempty(st.atom_constraints) ? collect(1:length(st.constraints)) :
            reduce(intersect, st.atom_constraints)
        some_atoms = length(whole) < length(st.constraints)
        constraints = st.constraints[whole]
        if some_atoms && isempty(whole)
            push!(results, (; structure=name, number=st.number_of_molecules, volume=NaN, error=NaN, density=NaN,
                note="region set by constraints of some atoms only"))
            continue
        end
        # The box to sample: the periodic cell, or the intersection of the bounding boxes
        # of the "inside" constraints
        lo = SVector{3,T}(-Inf, -Inf, -Inf)
        hi = SVector{3,T}(Inf, Inf, Inf)
        if has_pbc
            h = _cm_half_extent(system)
            lo, hi = center - h, center + h
        end
        for c in constraints
            b = _constraint_dispatch(_inside_box, c)
            isnothing(b) && continue
            lo, hi = max.(lo, b[1]), min.(hi, b[2])
        end
        note = some_atoms ? "upper bound: constraints of some atoms ignored" : ""
        if !all(isfinite, lo) || !all(isfinite, hi)
            push!(results, (; structure=name, number=st.number_of_molecules, volume=NaN, error=NaN, density=NaN,
                note="region not bounded (no periodic cell or inside constraint)" * (some_atoms ? "; constraints of some atoms ignored" : "")))
            continue
        end
        if any(hi .<= lo)
            push!(results, (; structure=name, number=st.number_of_molecules, volume=0.0, error=0.0, density=Inf,
                note="the constraints do not intersect"))
            continue
        end
        points = Vector{SVector{3,T}}(undef, npoints)
        extent = hi - lo
        for i in 1:npoints
            points[i] = lo + SVector{3,T}(rand(rng, T), rand(rng, T), rand(rng, T)) .* extent
        end
        # Points in the region: in the periodic cell (for non-orthorhombic cells, whose
        # bounding box is larger), and satisfying the constraints
        inside = trues(npoints)
        for i in 1:npoints
            x = points[i]
            if has_pbc && _wrap_cm(x, uc, center, system.wigner_seitz_cell) != x
                inside[i] = false
                continue
            end
            for c in constraints
                if _constraint_dispatch(constraint_penalty, c, x) > zero(T)
                    inside[i] = false
                    break
                end
            end
        end
        # Points within the tolerance of a fixed atom
        if !isempty(fixed)
            candidates = findall(inside)
            if !isempty(candidates)
                near = _near_fixed(points[candidates], fixed, tol, has_pbc, uc)
                for (k, i) in enumerate(candidates)
                    near[k] && (inside[i] = false)
                end
            end
        end
        p = count(inside) / npoints
        box_volume = prod(extent)
        volume = p * box_volume
        error = box_volume * sqrt(p * (1 - p) / npoints)
        density = st.number_of_molecules * _molar_mass(st.atoms) / _AVOGADRO / (volume * 1e-24)
        push!(results, (; structure=name, number=st.number_of_molecules, volume, error, density, note))
    end
    _print_free_volume(results, npoints)
    return results
end

# Whether each point is within `tol` of an atom of `fixed` (with the periodic images,
# if `has_pbc`)
#
# Without PBC, the cell lists are given a cell larger than the points and the fixed atoms
# by more than the cutoff, so that no periodic image is within it (as the packing does,
# see `_nonperiodic_unitcell`: the non-periodic mode of CellListMap does not compile with
# `juliac --trim`).
#
function _near_fixed(points::Vector{SVector{3,T}}, fixed::Vector{SVector{3,T}}, tol::T, has_pbc::Bool, uc) where {T}
    unitcell = if has_pbc
        SMatrix{3,3,T}(uc)
    else
        lo = reduce((a, b) -> min.(a, b), points; init=reduce((a, b) -> min.(a, b), fixed))
        hi = reduce((a, b) -> max.(a, b), points; init=reduce((a, b) -> max.(a, b), fixed))
        SMatrix{3,3,T}(Diagonal((hi - lo) .+ 3 * tol))
    end
    return _near_fixed_count(points, fixed, tol, unitcell) .> 0
end
function _near_fixed_count(points::Vector{SVector{3,T}}, fixed::Vector{SVector{3,T}}, tol::T, unitcell) where {T}
    sys = CellListMap.ParticleSystem(
        xpositions=points, ypositions=fixed, cutoff=tol, unitcell=unitcell,
        output=zeros(Int, length(points)), parallel=_nbatches() > 1, nbatches=_cl_nbatches(),
        validate_coordinates=_validate_positions,
    )
    pairwise!(sys) do pair, near
        pair.d < tol && (near[pair.i] += 1)
        return near
    end
    return sys.output
end

function _print_free_volume(results, npoints::Integer)
    _println()
    _println("  Volume available to each structure type (Monte Carlo, ", npoints, " points):")
    _println()
    _println("    ", _rpad("Structure", 24), _lpad("Number", 8), _lpad("Volume (Å³)", 16), _lpad("± error", 12), _lpad("Density (g/mL)", 16))
    for r in results
        if isnan(r.volume)
            _println("    ", _rpad(r.structure, 24), _lpad(string(r.number), 8), _lpad("-", 16), _lpad("-", 12), _lpad("-", 16), "   ", r.note)
        else
            _println("    ", _rpad(r.structure, 24), _lpad(string(r.number), 8), _lpad(_fmt(r.volume, 1), 16),
                _lpad(_fmt(r.error, 1), 12), _lpad(_fmt(r.density, 3), 16), isempty(r.note) ? "" : "   " * r.note)
        end
    end
    _println()
    _println("  The density is that at which the molecules of each structure type fill its region.")
    _println()
    return nothing
end

"""
    molecular_volume(atoms::AbstractVector{<:Atom}; probe = 1.4, resolution = 0.4)
    molecular_volume(input_file::String; probe = 1.4, resolution = 0.4)
    molecular_volume(system::PackmolSystem; probe = 1.4, resolution = 0.4)
    molecular_volume(molecules::AbstractVector{<:Pair}; probe = 1.4, resolution = 0.4)

The molecular volume (Å³): the volume from which a spherical probe of radius `probe` (Å)
is excluded by the van der Waals spheres of the atoms (the van der Waals volume and the
gaps between the atoms where the probe does not fit), computed on a grid of spacing
`resolution` (Å), as by MoloVol (Maglic & Lavendomme, J. Appl. Cryst. 55, 1033-1044 (2022), https://doi.org/10.1107/S1600576722004988).
Cavities where the probe fits are not part of the volume, whether they are connected to
the outside or not. With `probe = 0`, the van der Waals volume.

With the default probe (the radius of a water molecule), the molecular volume of a
protein is close to its volume in solution (that of its partial specific volume, about
0.73 mL/g). The radii are the van der Waals radii of the elements (`element_vdw_radius` of
PDBToolsBase), so that the volume is that of all-atom structures, whose elements are given
or follow from the names of the atoms. At the default resolution, the volume is within
about 1.5% of the converged one (slightly larger), at a fraction of a second for thousands
of atoms; the error decreases with the resolution, at a cost that grows as its inverse
cube. The probes are placed on the accessible surface and on its creases, where they
touch two or three atoms, which makes the volume converge fast with the resolution.

Of the atoms of a molecule, returns its volume. Of the structure types of an input file or
of a `PackmolSystem`, or of pairs of structure files and numbers of molecules (as
`["protein.pdb" => 1, "water.pdb" => 1000]`), returns `(; volume, molecules)`: the total
volume of the molecules, and for each structure type `(; structure, number, volume)`, the
volume of one molecule. Prints a summary.

# Example

```julia
using Packmol, PDBTools
molecular_volume(read_pdb("protein.pdb"))
molecular_volume(read_pdb("protein.pdb"); resolution = 0.2) # more precise, slower
```
"""
function molecular_volume(atoms::AbstractVector{<:Atom}; probe::Real=1.4, resolution::Real=0.4)
    probe >= 0 || throw(ArgumentError("\n\n    The probe radius must be non-negative, got $probe.\n\n"))
    resolution > 0 || throw(ArgumentError("\n\n    The resolution must be positive, got $resolution.\n\n"))
    isempty(atoms) && return 0.0
    x = SVector{3,Float64}[SVector{3,Float64}(a.x, a.y, a.z) for a in atoms]
    r = Float64[_vdw_radius(a) for a in atoms]
    return _molecular_volume(x, r, Float64(probe), Float64(resolution))
end

function molecular_volume(input_file::String; probe::Real=1.4, resolution::Real=0.4)
    return molecular_volume(read_packmol_input(input_file, Val(3), Float64); probe, resolution)
end

function molecular_volume(system::PackmolSystem; probe::Real=1.4, resolution::Real=0.4)
    names = String[]
    volumes = Float64[]
    numbers = Int[]
    for st in system.structure_types
        push!(names, basename(st.filename))
        push!(volumes, molecular_volume(st.atoms; probe, resolution))
        push!(numbers, st.number_of_molecules)
    end
    return _molecular_volume_summary(names, volumes, numbers, Float64(probe), Float64(resolution))
end

function molecular_volume(molecules::AbstractVector{<:Pair}; probe::Real=1.4, resolution::Real=0.4)
    names = String[]
    volumes = Float64[]
    numbers = Int[]
    for (molecule, number) in molecules
        file = expanduser(String(molecule))
        push!(names, basename(file))
        push!(volumes, molecular_volume(read_pdb(file); probe, resolution))
        push!(numbers, Int(number))
    end
    return _molecular_volume_summary(names, volumes, numbers, Float64(probe), Float64(resolution))
end

function _molecular_volume_summary(names::Vector{String}, volumes::Vector{Float64}, numbers::Vector{Int},
    probe::Float64, resolution::Float64)
    volume = sum(volumes[i] * numbers[i] for i in eachindex(volumes); init=0.0)
    _println()
    _println("  Molecular volumes (probe of ", _fmt(probe, 2), " Å, resolution of ", _fmt(resolution, 2), " Å):")
    _println()
    _println("    ", _rpad("Structure", 30), _lpad("Number", 9), _lpad("Volume (Å³)", 16), _lpad("Total (Å³)", 16))
    for i in eachindex(names)
        _println("    ", _rpad(names[i], 30), _lpad(string(numbers[i]), 9), _lpad(_fmt(volumes[i], 1), 16),
            _lpad(_fmt(volumes[i] * numbers[i], 1), 16))
    end
    _println()
    _println("  Total volume: ", _fmt(volume, 1), " Å³")
    _println()
    molecules = NamedTuple{(:structure, :number, :volume),Tuple{String,Int,Float64}}[
        (; structure=names[i], number=numbers[i], volume=volumes[i]) for i in eachindex(names)
    ]
    return (; volume, molecules)
end

# The van der Waals radius of the element of an atom
function _vdw_radius(atom::Atom)
    r = Float64(PDBToolsBase.element_vdw_radius(atom))
    if !(r > 0)
        throw(ArgumentError(string("\n\n    The van der Waals radius of atom ", atom.index, " (", atom.name,
            ") is unknown: set its element (columns 77-78 of PDB files).\n\n")))
    end
    return r
end

#
# The molecular volume of spheres at `x` of radii `r`, on a grid of spacing `g`: the grid
# points (voxel centers) that no point of a probe of radius `rp` reaches, with the center
# of the probe outside all spheres of radii `r + rp` (the probe-occupied space). The
# centers of the probe are on the surface of those spheres (the accessible surface): from
# each grid point outside the spheres and next to them, the nearest point of the surface,
# and the nearest points of its creases, where it touches two atoms (the circles where two
# spheres meet) or three (the cusps where three meet), which reach the bottoms of the
# grooves between the atoms. Without the creases, the probes miss slivers along them, of
# the width of the spacing: the volume converges linearly with the resolution (4-5% too
# large at 0.4 Å; 0.1-1.5% with them). With `rp = 0`, the van der Waals volume.
#
function _molecular_volume(x::Vector{SVector{3,Float64}}, r::Vector{Float64}, rp::Float64, g::Float64)
    R = r .+ rp
    pad = maximum(R) + rp + 2g
    lo = reduce((a, b) -> min.(a, b), x) .- pad
    hi = reduce((a, b) -> max.(a, b), x) .+ pad
    n = ceil.(Int, (hi - lo) ./ g)
    # The space excluded to the centers of the probe
    core = falses(n[1], n[2], n[3])
    _stamp_spheres!(core, x, R, lo, g)
    rp == 0 && return count(core) * g^3
    # The centers of the probe next to that space, moved to its surface and creases
    centers = SVector{3,Float64}[]
    @inbounds for k in 2:n[3]-1, j in 2:n[2]-1, i in 2:n[1]-1
        core[i, j, k] && continue
        if core[i-1, j, k] || core[i+1, j, k] || core[i, j-1, k] || core[i, j+1, k] || core[i, j, k-1] || core[i, j, k+1]
            push!(centers, lo + (SVector{3,Float64}(i, j, k) .- 0.5) .* g)
        end
    end
    _probe_centers_on_surface!(centers, x, R, lo, g)
    # The space reached by the probe: outside the excluded space, or within the probes
    reached = .!core
    _stamp_spheres!(reached, centers, fill(rp, length(centers)), lo, g)
    return count(!, reached) * g^3
end

# Mark the grid points (voxel centers) within spheres at `x` of radii `r`
function _stamp_spheres!(grid::BitArray{3}, x::Vector{SVector{3,Float64}}, r::Vector{Float64},
    lo::SVector{3,Float64}, g::Float64)
    n = size(grid)
    @inbounds for (xi, ri) in zip(x, r)
        imin = max.(1, floor.(Int, (xi .- ri .- lo) ./ g))
        imax = min.(SVector(n), ceil.(Int, (xi .+ ri .- lo) ./ g) .+ 1)
        r2 = ri^2
        for k in imin[3]:imax[3], j in imin[2]:imax[2], i in imin[1]:imax[1]
            p = lo + (SVector{3,Float64}(i, j, k) .- 0.5) .* g
            sum(abs2, p - xi) <= r2 && (grid[i, j, k] = true)
        end
    end
    return grid
end

# Move each center of the probe (a grid point outside the spheres `R`, next to them) to
# the nearest point of the surface of the spheres (that of the nearest sphere, which is
# outside all others), and add the nearest points of the creases of the surface: the
# point of the circle where the nearest two spheres meet, and the point where the
# nearest three meet, if they are near (within 2g) and outside all other spheres
function _probe_centers_on_surface!(centers::Vector{SVector{3,Float64}}, x::Vector{SVector{3,Float64}},
    R::Vector{Float64}, lo::SVector{3,Float64}, g::Float64)
    # Cells of the atoms, of side larger than the radii of the spheres and the distances
    # of the centers to them: the spheres near a center are those of the cells within two
    side = maximum(R) + 2g
    nc = SVector{3,Int}(ceil.(Int, (reduce((a, b) -> max.(a, b), x) - lo) ./ side) .+ 3)
    cells = [Int[] for _ in 1:prod(nc)]
    lin = LinearIndices((nc[1], nc[2], nc[3]))
    for (i, xi) in enumerate(x)
        c = floor.(Int, (xi - lo) ./ side) .+ 1
        push!(cells[lin[c[1], c[2], c[3]]], i)
    end
    near = Int[]
    creases = SVector{3,Float64}[]
    for ic in eachindex(centers)
        p = centers[ic]
        c = floor.(Int, (p - lo) ./ side) .+ 1
        empty!(near)
        for dk in -2:2, dj in -2:2, di in -2:2
            a, b, d = c[1] + di, c[2] + dj, c[3] + dk
            (1 <= a <= nc[1] && 1 <= b <= nc[2] && 1 <= d <= nc[3]) || continue
            append!(near, cells[lin[a, b, d]])
        end
        isempty(near) && continue
        # The three nearest spheres (distances to their surfaces)
        i1, i2, i3 = 0, 0, 0
        d1, d2, d3 = Inf, Inf, Inf
        for i in near
            dist = norm(p - x[i]) - R[i]
            if dist < d1
                i3, d3 = i2, d2
                i2, d2 = i1, d1
                i1, d1 = i, dist
            elseif dist < d2
                i3, d3 = i2, d2
                i2, d2 = i, dist
            elseif dist < d3
                i3, d3 = i, dist
            end
        end
        centers[ic] = x[i1] + R[i1] * (p - x[i1]) / norm(p - x[i1])
        i2 == 0 && continue
        found, q = _crease_point(p, x[i1], R[i1], x[i2], R[i2])
        if found && norm(q - p) < 2g && !_inside_spheres(q, x, R, near, i1, i2, 0)
            push!(creases, q)
        end
        i3 == 0 && continue
        found, q = _cusp_point(p, x[i1], R[i1], x[i2], R[i2], x[i3], R[i3])
        if found && norm(q - p) < 2g && !_inside_spheres(q, x, R, near, i1, i2, i3)
            push!(creases, q)
        end
    end
    append!(centers, creases)
    return centers
end

# Whether `q` is inside one of the spheres `near`, other than `i`, `j` and `k`
function _inside_spheres(q::SVector{3,Float64}, x::Vector{SVector{3,Float64}}, R::Vector{Float64},
    near::Vector{Int}, i::Int, j::Int, k::Int)
    for l in near
        (l == i || l == j || l == k) && continue
        sum(abs2, q - x[l]) < (R[l] - 1e-9)^2 && return true
    end
    return false
end

# The point nearest to `p` of the circle where the spheres (xi, Ri) and (xj, Rj) meet
# (`false` as the first value if they do not meet)
function _crease_point(p::SVector{3,Float64}, xi::SVector{3,Float64}, Ri::Float64,
    xj::SVector{3,Float64}, Rj::Float64)
    v = xj - xi
    d = norm(v)
    (d >= Ri + Rj || d <= abs(Ri - Rj)) && return false, p
    n = v / d
    a = (d^2 + Ri^2 - Rj^2) / (2d)
    h = sqrt(max(Ri^2 - a^2, 0.0))
    center = xi + a * n
    w = p - center
    w = w - dot(w, n) * n
    nw = norm(w)
    nw < 1e-12 && return false, p
    return true, center + h * w / nw
end

# The point nearest to `p` of the two where the spheres (x1, R1), (x2, R2) and (x3, R3)
# meet (`false` as the first value if they do not meet)
function _cusp_point(p::SVector{3,Float64}, x1::SVector{3,Float64}, R1::Float64,
    x2::SVector{3,Float64}, R2::Float64, x3::SVector{3,Float64}, R3::Float64)
    ex = x2 - x1
    d = norm(ex)
    d < 1e-12 && return false, p
    ex = ex / d
    i = dot(ex, x3 - x1)
    ey = x3 - x1 - i * ex
    ney = norm(ey)
    ney < 1e-12 && return false, p
    ey = ey / ney
    ez = cross(ex, ey)
    j = dot(ey, x3 - x1)
    X = (R1^2 - R2^2 + d^2) / (2d)
    Y = (R1^2 - R3^2 + i^2 + j^2) / (2j) - i / j * X
    z2 = R1^2 - X^2 - Y^2
    z2 < 0 && return false, p
    z = sqrt(z2)
    base = x1 + X * ex + Y * ey
    a = base + z * ez
    b = base - z * ez
    return true, (norm(a - p) < norm(b - p) ? a : b)
end
