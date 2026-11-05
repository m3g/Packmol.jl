#
# Tools for building systems by hand (public, not exported): the cross-section of a
# molecule in a slab, and the data of the built-in lipids.
#
public cross_section, builtin_lipid

"""
    Packmol.cross_section(atoms; z = (-Inf, Inf), radius = 2.0, spacing = 0.25)

The area (Å²) of the projection on the `xy` plane of the atoms with `z[1] ≤ z ≤ z[2]`:
the union of the discs of radius `radius` around them, computed on a grid of spacing
`spacing`. `atoms` is a vector of `Atom`s (as read by `PDBTools.read_pdb`) or of
coordinates.

It is the area that a membrane protein takes from the lipids of a leaflet: with the
protein oriented along `z` and the middle of the bilayer at `z = 0` (as the structures
of the [OPM database](https://opm.phar.umich.edu)), the number of lipids of the upper
leaflet of a membrane of area `A` is `(A - cross_section(protein; z = (0, h))) /
area_per_lipid`, with `h` the half-thickness of the membrane (about 15 to 17 Å for the
hydrophobic core). With `radius` equal to the tolerance of the packing, it is the area
that the atoms of the lipids cannot occupy.

# Example

```julia
using Packmol, PDBTools
protein = read_pdb("protein.pdb")      # oriented as in OPM
Packmol.cross_section(protein; z = (0, 17))
```
"""
function cross_section(atoms::AbstractVector; z = (-Inf, Inf), radius::Real = 2.0, spacing::Real = 0.25)
    xs = [_coordinates(a) for a in atoms]
    xs = [x for x in xs if z[1] <= x[3] <= z[2]]
    isempty(xs) && return 0.0
    lo = reduce((a, b) -> min.(a, b), xs) .- radius .- spacing
    hi = reduce((a, b) -> max.(a, b), xs) .+ radius .+ spacing
    return _cross_section(xs, lo, hi, Float64(radius); h = Float64(spacing))
end
_coordinates(a::Atom) = SVector{3,Float64}(a.x, a.y, a.z)
_coordinates(x::AbstractVector) = SVector{3,Float64}(x[1], x[2], x[3])

"""
    Packmol.builtin_lipid(name)

The data of a built-in lipid of the [Membrane recipe](@ref "Membrane system"), given by
its name (`"POPC"`, `"POPE"`, `"POPS"`, `"POPG"`, `"DPPC"`, `"DOPC"` or `"DMPC"`, case
insensitive), as a named tuple:

- `file`: the structure file (PDB, CHARMM36 atom names);
- `head`: the index of the atom of the polar head (the nitrogen, or the phosphorus);
- `chain_ends`: the indices of the last carbon atoms of the chains;
- `area`: the area per lipid of the actual membrane (Å²);
- `length`: the distance from the head atom to the middle of the chain ends (Å).

These are the data needed to build membranes of other shapes (vesicles, curved or
asymmetric bilayers) with constraints on the heads and on the chain ends: in a flat
bilayer, the `Membrane` recipe places the chain ends within 6 Å of the middle, and the
heads beyond `6 + 0.85 * length` (see [Membrane system](@ref)).

# Example

```julia-repl
julia> lipid = Packmol.builtin_lipid("POPC");

julia> lipid.head, lipid.chain_ends, lipid.area
(1, [88, 131], 64.0)
```
"""
function builtin_lipid(name::AbstractString)
    lipid = _builtin_lipid(name)
    isnothing(lipid) && throw(ArgumentError("""\n
        $name is not a built-in lipid. The built-in lipids are: $(join(sort(collect(keys(_BUILTIN_LIPIDS))), ", ")).

    """))
    file = joinpath(_MEMBRANE_LIPIDS_DIR, lowercase(name) * ".pdb")
    x = [_coordinates(a) for a in read_pdb(file)]
    ends = sum(x[lipid.chain_ends]) / length(lipid.chain_ends)
    return (file = file, head = lipid.head, chain_ends = copy(lipid.chain_ends), area = lipid.area,
        length = norm(x[lipid.head] - ends))
end
