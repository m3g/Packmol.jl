# Restarts (`restart_from`) use coordinate files (PDB, XYZ or mmCIF, by the file
# extension): the position of each molecule is recovered from the coordinates (see
# `_align_molecule`), of the output of a previous run or of a simulation. The restart
# files of Fortran Packmol (raw center of mass and Euler angles) are not supported.
#
_restart_file_format(path::AbstractString) = _file_format(path, "")

function _check_restart_file(path::AbstractString, keyword::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext in (".pdb", ".xyz", ".cif", ".mmcif") && return path
    throw(ArgumentError("""\n
        $keyword $path: restart files must be coordinate files (.pdb, .xyz or .cif).
        The restart files of Fortran Packmol are not supported: restart from the
        coordinates of a previous run (its output file).

    """))
end
_check_restart_file(::Nothing, ::AbstractString) = nothing
_check_restart_file(::AbstractVector{<:Atom}, ::AbstractString) = nothing

function _check_restart_files(packmol_system::PackmolSystem)
    _check_restart_file(packmol_system.restart_from, "restart_from")
    for st in packmol_system.structure_types
        _check_restart_file(st.restart_from, "restart_from")
    end
    return nothing
end

#
# Recovers a molecule's (cm, angles) from its atoms' observed Cartesian
# positions, by rigid-body-aligning `reference` (that structure type's own
# `reference_coordinates` — already centered at its own centroid) onto
# `observed` (raw, uncentered) via the Kabsch algorithm: `cm` is `observed`'s
# centroid, and `angles` is the least-squares-optimal rotation taking
# `reference` onto `observed - cm`, decomposed back into `eulermat`'s
# parameterization.
#
function _align_molecule(reference::AbstractVector{SVector{3,T}}, observed::AbstractVector{SVector{3,T}}) where {T}
    cm = sum(observed) / length(observed)
    R = kabsch_rotation(reference, observed .- Ref(cm))
    return MoleculePosition(cm, euler_angles(R))
end

#
# Reads restart positions for a sequence of structure-type "segments" — each
# a (number_of_molecules, natoms, reference_coordinates, extract) tuple, in
# the same order the atoms appear in `atoms` — by rigid-body-aligning every
# molecule's atoms against its own type's template (see `_align_molecule`).
# Used for both a coordinate-file restart_from and an in-memory `Vector{<:Atom}`
# one; both go through this same function once the atoms are in hand.
#
# `extract` lets a segment's atoms be *skipped* (counted toward the expected
# total and advanced over, but not aligned/returned) rather than always
# producing a position: a whole-system restart_from PDB is whatever
# `write_output` last wrote for the *entire* system, fixed molecules
# included (interleaved with free ones in structure-declaration order) — so
# the expected atom count and the file's actual atom-by-atom layout both
# need every structure type accounted for, even though a fixed molecule's
# own position is never read back (it's already pinned by its own `fixed`
# keyword regardless of what the restart source says).
#
function _restart_positions_from_atoms(
    atoms::AbstractVector{<:Atom}, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}},Bool}},
) where {D,T}
    expected = sum(nmols * natoms for (nmols, natoms, _, _) in segments)
    length(atoms) == expected || error(
        "Restart coordinates have $(length(atoms)) atoms, expected $expected " *
        "($(join(("$nmols × $natoms" for (nmols, natoms, _, _) in segments), " + ")))."
    )
    positions = MoleculePosition{D,T}[]
    iat = 0
    for (nmols, natoms, refcoords, extract) in segments
        for _ in 1:nmols
            if extract
                observed = SVector{D,T}[_atom_position(atoms[iat+j], Val(D), T) for j in 1:natoms]
                push!(positions, _align_molecule(refcoords, observed))
            end
            iat += natoms
        end
    end
    return positions
end

#
# Restart positions for one or more structure-type `segments` (see
# `_restart_positions_from_atoms` for the tuple layout and what `extract`
# means), from an in-memory atom vector (Julia API only) or a coordinate file.
#
# A restart source given as atoms (e.g. the `Structure` returned by `read_pdb`) is
# stored as a `Vector{Atom{Nothing}}`.
_restart_source(source) = source
_restart_source(source::AbstractVector{<:Atom}) = Vector{Atom{Nothing}}(source)

# How a restart source is reported in the output.
_restart_source_name(source::AbstractString) = source
_restart_source_name(source::AbstractVector{<:Atom}) = "atom vector ($(length(source)) atoms)"

function _restart_positions(source::AbstractVector{<:Atom}, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}},Bool}}, ::Type{MoleculePosition{D,T}}) where {D,T}
    _restart_positions_from_atoms(source, segments)
end
function _restart_positions(source::AbstractString, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}},Bool}}, ::Type{MoleculePosition{D,T}}) where {D,T}
    atoms, _ = _read_structure_file(source, _restart_file_format(source))
    _restart_positions_from_atoms(atoms, segments)
end

#
