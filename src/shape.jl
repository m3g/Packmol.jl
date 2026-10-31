#
# Shape of molecules: descriptors from the gyration tensor, to choose the conformations
# of a molecule that pack best (rod-like for membranes, sphere-like for homogeneous
# mixtures).
#

"""
    shape_descriptors(atoms_or_coordinates)
    shape_descriptors(filename)

Descriptors of the shape of a molecule, from the eigenvalues `λ₁ ≥ λ₂ ≥ λ₃` of its
gyration tensor (the variances of the coordinates of its atoms along its principal
axes, all atoms with the same weight). Returns a named tuple with:

- `principal_variances`: `(λ₁, λ₂, λ₃)`, in Å²;
- `radius_of_gyration`: `√(λ₁ + λ₂ + λ₃)`, in Å;
- `asphericity`: `λ₁ - (λ₂ + λ₃)/2`, in Å²: how much the first principal axis
  dominates the others. It is zero for a spherically symmetric molecule and greatest,
  among conformations of a molecule, for the most elongated (rod-like) one;
- `acylindricity`: `λ₂ - λ₃`, in Å²: zero for a molecule with cylindrical symmetry;
- `anisotropy`: the relative shape anisotropy, `(b² + 3c²/4)/R_g⁴` (`b` the asphericity,
  `c` the acylindricity), from 0 (sphere-like) to 1 (atoms on a line); 1/4 for a flat
  disc.

The argument is a vector of atoms or of coordinates, or a structure file. For a file
with several models (conformations, see [`structure_type`](@ref)), a vector with the
descriptors of each is returned.

Molecules are rigid in Packmol: the shape of the conformations affects how densely they
pack. Rod-like conformations (large asphericity) pack best in membranes, and
sphere-like conformations (small anisotropy) in homogeneous mixtures. See
[`select_conformations`](@ref).

# Example

```julia-repl
julia> shape_descriptors("water.pdb").anisotropy
```
"""
function shape_descriptors(x::AbstractVector{<:AbstractVector{<:Real}})
    isempty(x) && throw(ArgumentError("shape_descriptors: no coordinates given."))
    D = length(first(x))
    c = sum(SVector{3,Float64}(v[1], v[2], D == 3 ? v[3] : 0.0) for v in x) / length(x)
    S = zeros(SMatrix{3,3,Float64})
    for v in x
        r = SVector{3,Float64}(v[1], v[2], D == 3 ? v[3] : 0.0) - c
        S += r * r'
    end
    S /= length(x)
    λ = sort(max.(eigvals(Symmetric(Matrix(S))), 0.0); rev=true)
    rg2 = sum(λ)
    b = λ[1] - (λ[2] + λ[3]) / 2
    c2 = λ[2] - λ[3]
    anisotropy = rg2 > 0 ? (b^2 + 3c2^2 / 4) / rg2^2 : 0.0
    return (
        principal_variances=(λ[1], λ[2], λ[3]),
        radius_of_gyration=sqrt(rg2),
        asphericity=b,
        acylindricity=c2,
        anisotropy=anisotropy,
    )
end
shape_descriptors(atoms::AbstractVector{<:Atom}) =
    shape_descriptors([SVector{3,Float64}(at.x, at.y, at.z) for at in atoms])
function shape_descriptors(filename::String; filetype::String="pdb")
    conformations, _, _ = _read_conformations(filename, filetype)
    length(conformations) == 1 && return shape_descriptors(only(conformations))
    return [shape_descriptors(conf) for conf in conformations]
end

"""
    select_conformations(input, output; n=1, shape=:rod, filetype="pdb")

Write to `output` a structure file with the `n` (by default, one) conformations of a
molecule, from those of `input` (a structure file with several models, or a vector of
files), that are most rod-like (`shape = :rod`: greatest asphericity) or most
sphere-like (`shape = :sphere`: smallest relative shape anisotropy), as models, the
best first. The output file can be
used in a structure block of an input file, or in [`structure_type`](@ref), for a
molecule packed in these conformations. Returns the indices of the selected
conformations in `input`, and their shape descriptors (see [`shape_descriptors`](@ref)).

The conformations can be the frames of a simulation of a single molecule (in vacuum, at
a high temperature, to sample many conformations), saved as the models of a PDB file,
or the molecules of an equilibrated system.

Rod-like conformations pack best in membranes: a POPC bilayer packs at its actual area
per lipid (64 Å²) from copies of the most rod-like conformation of an equilibrated
membrane, but not from all its conformations. Sphere-like conformations pack best in
homogeneous mixtures. A single conformation packs best: the same bilayer packs less
often from the 5 most rod-like conformations, and even less from the 10 most rod-like.

# Example

```julia-repl
julia> index, descriptors = select_conformations("popc_conformations.pdb", "popc_rod.pdb")
```
"""
function select_conformations(
    input::Union{String,AbstractVector{String}}, output::String;
    n::Int=1, shape::Symbol=:rod, filetype::String="pdb",
)
    shape in (:rod, :sphere) || throw(ArgumentError("""\n
        select_conformations: shape must be :rod or :sphere, got :$shape.

    """))
    conformations, _, _ = _read_conformations(input, filetype)
    n >= 1 || throw(ArgumentError("select_conformations: n must be at least 1."))
    n = min(n, length(conformations))
    descriptors = [shape_descriptors(conf) for conf in conformations]
    order = shape == :rod ?
        sortperm(descriptors; by=d -> d.asphericity, rev=true) :
        sortperm(descriptors; by=d -> d.anisotropy)
    selected = order[1:n]
    open(output, "w") do io
        println(io, "REMARK   $n conformations selected by select_conformations (shape = $shape)")
    end
    for (k, i) in enumerate(selected)
        write_pdb(output, conformations[i];
            append=true, header="MODEL     " * lpad(k, 4), footer="ENDMDL", connectivity=false,
        )
    end
    open(io -> println(io, "END"), output, "a")
    return selected, descriptors[selected]
end
