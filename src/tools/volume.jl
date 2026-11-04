#
# Two estimates for setting up a system: the volume that a set of molecules fills at a
# given density (`guess_volume`), and the volume of the region available to the
# molecules of each structure type of an input, given its constraints and the fixed
# molecules (`free_volume`, by Monte Carlo).
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

#
# The bounding box of an "inside" constraint (`nothing` for the others)
#
_inside_box(c::Box{Inside}) = (c.center - c.sides / 2, c.center + c.sides / 2)
_inside_box(c::Cube{Inside,T}) where {T} = (c.center .- c.side / 2, c.center .+ c.side / 2)
_inside_box(c::Sphere{Inside}) = (c.center .- c.radius, c.center .+ c.radius)
function _inside_box(c::Cylinder{Inside})
    tip = c.center + c.axis / norm(c.axis) * c.length
    return (min.(c.center, tip) .- c.radius, max.(c.center, tip) .+ c.radius)
end
function _inside_box(c::Ellipsoid{Inside})
    half = abs(c.scale) * SVector(c.a, c.b, c.c)
    return (c.center - half, c.center + half)
end
_inside_box(::Constraint) = nothing

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
