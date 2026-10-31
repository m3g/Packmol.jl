#
# Residue number for one atom, following Fortran Packmol's `resnumbers` modes
# (`output.f90`): 0 renumbers each molecule sequentially (constant 1 for a
# fixed molecule, since there's only one to number), 1 keeps the atom's
# residue number from its own template PDB verbatim, 2 offsets the template's
# residue numbers by a running counter so consecutive molecules of the same
# type get consecutive residues (matching multi-residue templates, e.g. a
# solvated chain of amino acids), and 3 gives every molecule in the whole
# system (across all structure types) its own residue number. All modes wrap
# at 9999 (the PDB residue-number field's width) with 0 mapping to 9999, to
# match the Fortran output exactly.
#
function _resnum(mode::Int, is_fixed::Bool, imol_local::Int, imark::Integer, ifres::Integer, irescount::Int, iimol::Int)
    r = if mode == 0
        is_fixed ? 1 : mod(imol_local, 9999)
    elseif mode == 1
        imark
    elseif mode == 2
        mod(imark - ifres + irescount, 9999)
    elseif mode == 3
        mod(iimol, 9999)
    else
        error("Unknown resnumbers mode $mode (must be 0, 1, 2, or 3).")
    end
    return r == 0 ? 9999 : r
end

#
# Chain letter for auto-assigned (non-fixed `chain` keyword) chains, following
# Fortran Packmol's `chainc` subroutine: 1-26 -> 'A'-'Z', 27-35 -> '1'-'9',
# 36 -> '0', beyond that there are no more identifiers left so it falls back
# to '#' (Fortran's own sentinel for "ran out").
#
function _chainc(i::Int)
    i == 0 && return ' '
    1 <= i <= 26 && return Char(UInt8('A') + i - 1)
    27 <= i <= 35 && return Char(UInt8('1') + i - 27)
    i == 36 && return '0'
    return '#'
end

#
# (a, b, c, α, β, γ) in the CRYST1/PDB convention from a 3-column unitcell
# matrix (columns are the cell vectors, as built by `unitcell_matrix`): a, b,
# c are the column norms; α is the angle between b and c, β between a and c,
# γ between a and b (all in degrees).
#
function _unitcell_abc_angles(unitcell::AbstractMatrix{T}) where {T}
    D = size(unitcell, 1)
    a_vec = unitcell[:, 1]
    b_vec = unitcell[:, 2]
    a = norm(a_vec)
    b = norm(b_vec)
    if D == 2
        c = one(T)
        α = β = T(90)
        γ = acosd(clamp(dot(a_vec, b_vec) / (a * b), -one(T), one(T)))
    else
        c_vec = unitcell[:, 3]
        c = norm(c_vec)
        α = acosd(clamp(dot(b_vec, c_vec) / (b * c), -one(T), one(T)))
        β = acosd(clamp(dot(a_vec, c_vec) / (a * c), -one(T), one(T)))
        γ = acosd(clamp(dot(a_vec, b_vec) / (a * b), -one(T), one(T)))
    end
    return a, b, c, α, β, γ
end

#
# CRYST1 (+ explanatory REMARK, when the box isn't an explicit PBC unitcell)
# header lines, or an empty vector if neither `add_box_sides` nor PBC is set
# (matching Fortran's `if(add_box_sides .or. using_pbc)` gate). Precision is
# %9.3f/%7.2f (the standard PDB CRYST1 field widths) rather than Fortran
# Packmol's own %9.2f/%7.2f — a strictly more precise, still spec-compliant
# choice, not a behavioral difference worth preserving.
#
function _cryst1_lines(packmol_system::PackmolSystem{D,T}, atom_positions) where {D,T}
    add_box_sides = packmol_system.add_box_sides
    has_pbc = !isnothing(packmol_system.unitcell)
    (!add_box_sides && !has_pbc) && return String[]

    lines = String[]
    if has_pbc
        a, b, c, α, β, γ = _unitcell_abc_angles(packmol_system.unitcell)
    else
        push!(lines, "REMARK  CRYST1 info below is (extrema(coordinates) +/- 1.1*tolerance) because no explicit")
        push!(lines, "REMARK  PBCs were defined. To apply PBCs, use the `pbc` keyword.")
        lo, hi = compute_bounding_box(atom_positions)
        margin = T(1.1) * packmol_system.tolerance
        sides = (hi .- lo) .+ 2margin
        a, b = sides[1], sides[2]
        c = D == 3 ? sides[3] : one(T)
        α = β = γ = T(90)
    end
    push!(lines, @sprintf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d", a, b, c, α, β, γ, "P 1", 1))
    return lines
end

#
# Restarts (`restart_from`/`restart_to`) use coordinate files (PDB, XYZ or mmCIF, by
# the file extension): `restart_to` writes the coordinates of the system (or of one
# structure type) every time the output is written, and `restart_from` recovers the
# position of each molecule from coordinates (see `_align_molecule`). The restart
# files of Fortran Packmol (raw center of mass and Euler angles) are not supported.
#
_restart_file_format(path::AbstractString) = _file_format(path, "")

function _check_restart_file(path::AbstractString, keyword::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext in (".pdb", ".xyz", ".cif", ".mmcif") && return path
    throw(ArgumentError("""\n
        $keyword $path: restart files must be coordinate files (.pdb, .xyz or .cif).
        The restart files of Fortran Packmol are not supported: restart from the
        coordinates of a previous run (its output file, or a restart_to file).

    """))
end
_check_restart_file(::Nothing, ::AbstractString) = nothing
_check_restart_file(::AbstractVector{<:Atom}, ::AbstractString) = nothing

function _check_restart_files(packmol_system::PackmolSystem)
    _check_restart_file(packmol_system.restart_from, "restart_from")
    _check_restart_file(packmol_system.restart_to, "restart_to")
    for st in packmol_system.structure_types
        _check_restart_file(st.restart_from, "restart_from")
        _check_restart_file(st.restart_to, "restart_to")
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
# Write any configured restart_to files (whole-system and/or per structure
# type): the coordinates of `output_atoms` (the atoms of `get_atoms`), or of
# the atoms of one structure type.
#
function _write_restart_files(packmol_system::PackmolSystem, output_atoms::Vector{Atom{Nothing}})
    restart_to = packmol_system.restart_to
    isnothing(restart_to) || _write_coordinates(restart_to, output_atoms)
    iat_offset = 0
    for st in packmol_system.structure_types
        natoms = _ncopies(st) * st.natoms
        st_restart_to = st.restart_to
        if !isnothing(st_restart_to)
            _write_coordinates(st_restart_to, output_atoms[iat_offset+1:iat_offset+natoms])
        end
        iat_offset += natoms
    end
    return nothing
end

function _write_coordinates(filename::String, atoms::Vector{Atom{Nothing}})
    if _restart_file_format(filename) == "xyz"
        write_xyz(filename, Structure(atoms; title="Built with Packmol"))
    else
        write_pdb(filename, atoms)
    end
    return filename
end

#
# Per-molecule rigid CM shift applied to `atom_positions` in place, wrapping
# each molecule's center of mass into the unit cell (centered at
# `unitcell_center`) as a single rigid unit — every atom of that molecule
# carried along by the same offset — rather than wrapping each atom's
# absolute position independently. Wrapping atoms independently lets a
# molecule straddling a periodic boundary be torn in two (part of it wrapped
# to the opposite face while the rest stays put), which can land those
# wrapped atoms on top of whatever else sits there — visible as overlapping
# atoms in the output even though the packing itself converged. See the
# analogous comment on `_constraint_fg!` in interatomic_distance_fg.jl for
# the same fix applied during optimization.
#
# Which fundamental-domain shape to wrap into is `wrap`: `wrap_to_center`
# (the ordinary skewed-parallelepiped shape) or `triclinic_to_dodecahedral`
# (the rhombic-dodecahedron/Wigner-Seitz "compact" shape) — chosen by
# `get_atoms` from `packmol_system.periodic_boundary_style`.
#
function _wrap_molecules!(
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
    unitcell::AbstractMatrix,
    center::SVector{D,T},
    wrap::F,
) where {D,T,F}
    imol_offset = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        natoms_st = st.natoms
        nmols_st = st.number_of_molecules
        # Fixed molecules are written where they were placed (they are inside the cell)
        if st.fixed.fixed
            imol_offset += nmols_st
            iat_offset += nmols_st * natoms_st
            continue
        end
        for i in 1:nmols_st
            imol = imol_offset + i
            mp = packmol_system.molecule_positions[imol]
            shift = wrap(mp.cm, unitcell, center) - mp.cm
            if !iszero(shift)
                iat_first = iat_offset + (i - 1) * natoms_st
                for j in 1:natoms_st
                    atom_positions[iat_first + j] += shift
                end
            end
        end
        imol_offset += nmols_st
        iat_offset += nmols_st * natoms_st
    end
    return atom_positions
end

"""
    get_atoms(packmol_system::PackmolSystem)

Returns the atoms of the packed system, with their current coordinates.

The result is a `Vector{PDBToolsBase.Atom}`, generated from
`packmol_system.molecule_positions` on each call (never cached on `packmol_system`
itself), the same way `write_output` computes what it writes to file. Residue numbering follows each structure
type's `residue_numbering` (`resnumbers`) setting, and chain identifiers
follow `chain`/`changechains`.

When periodic boundary conditions are set (`packmol_system.unitcell`), each
molecule's center of mass is wrapped into the unit cell as a rigid unit (a
molecule straddling a periodic boundary is never torn in two across it), in
the representation selected by `packmol_system.periodic_boundary_style`:
the ordinary skewed-parallelepiped shape (`:triclinic`, the default), or,
after [`triclinic_to_dodecahedral`](@ref)/[`triclinic_to_octahedral`](@ref)
`(packmol_system)` has been called, the corresponding Wigner-Seitz "compact"
shape (`:dodecahedral`/`:octahedral`) — see those functions and their
`dodecahedral_to_triclinic`/`octahedral_to_triclinic` counterparts.
"""
function get_atoms(packmol_system::PackmolSystem{D,T}) where {D,T}
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    unitcell, unitcell_center = packmol_system.unitcell, packmol_system.unitcell_center
    if !isnothing(unitcell) && !isnothing(unitcell_center)
        # `:dodecahedral`/`:octahedral` share the very same "nearest
        # periodic image" search (see periodic_cells.jl) — only the
        # `unitcell` matrix itself differs between the two box shapes.
        if packmol_system.periodic_boundary_style === :triclinic
            _wrap_molecules!(atom_positions, packmol_system, unitcell, unitcell_center, wrap_to_center)
        else
            _wrap_molecules!(atom_positions, packmol_system, unitcell, unitcell_center, _nearest_periodic_image)
        end
    end

    output_atoms = Atom{Nothing}[]
    iat = 0
    atom_serial = 0
    iimol = 0       # molecule index across the whole system (fixed + free) — resnumbers mode 3
    irescount = 1   # running residue-number offset across the whole system — resnumbers mode 2
    ichain = 0      # count of auto-assigned chain letter groups so far, across all structure types
    # The conformations of a molecule (consecutive structure types, see `conformation`)
    # are numbered, and get chains, as a single structure type: `imol_offset` is the
    # number of molecules of the previous conformations of the set
    imol_offset = 0
    odd_chain = even_chain = ' '
    for st in packmol_system.structure_types
        st.conformation[1] == 1 && (imol_offset = 0)
        ifres = st.atoms[1].resnum
        ilres = st.atoms[end].resnum
        nres = ilres - ifres + 1
        resnumbers_mode = st.residue_numbering == -1 ? (nres == 1 ? 0 : 1) : st.residue_numbering

        if st.fixed.fixed
            iimol += 1
            for j in 1:st.natoms
                iat += 1
                atom_serial += 1
                atom = copy(st.atoms[j])
                pos = atom_positions[iat]
                atom.index = atom_serial
                atom.x = pos[1]
                atom.y = pos[2]
                atom.z = D == 3 ? pos[3] : zero(eltype(pos))
                atom.resnum = _resnum(resnumbers_mode, true, 1, st.atoms[j].resnum, ifres, irescount, iimol)
                # A fixed molecule with no explicit `chain` keeps its own
                # template chain id (there's no "auto-cycle per molecule" to
                # do for a single fixed copy).
                atom.chain = isnothing(st.chain) ? st.atoms[j].chain : string(st.chain)
                _set_segid!(atom, st.segid)
                push!(output_atoms, atom)
            end
            irescount += nres
        else
            for imol_in_type in 1:st.number_of_molecules
                imol_local = imol_offset + imol_in_type
                iimol += 1
                if isnothing(st.chain)
                    if imol_local == 1 || mod(imol_local, 9999) == 1
                        ichain += 1
                        odd_chain = _chainc(ichain)
                        if st.changechains
                            ichain += 1
                            even_chain = _chainc(ichain)
                        else
                            even_chain = odd_chain
                        end
                    end
                    write_chain = iseven(imol_local) ? even_chain : odd_chain
                else
                    write_chain = st.chain
                end
                for j in 1:st.natoms
                    iat += 1
                    atom_serial += 1
                    atom = copy(st.atoms[j])
                    pos = atom_positions[iat]
                    atom.index = atom_serial
                    atom.x = pos[1]
                    atom.y = pos[2]
                    atom.z = D == 3 ? pos[3] : zero(eltype(pos))
                    atom.resnum = _resnum(resnumbers_mode, false, imol_local, st.atoms[j].resnum, ifres, irescount, iimol)
                    atom.chain = string(write_chain)
                    _set_segid!(atom, st.segid)
                    push!(output_atoms, atom)
                end
                irescount += nres
            end
            imol_offset += st.number_of_molecules
        end
    end

    return output_atoms
end

_set_segid!(::Atom, ::Nothing) = nothing
_set_segid!(atom::Atom, segid::String) = (atom.segname = segid; nothing)

#
# The same positions `get_atoms` just computed, repackaged as `SVector{D,T}`
# from the `Atom`s it returned — cheap (O(natoms)), and avoids computing
# atom positions a second time just for `_cryst1_lines`'s no-PBC bounding-box
# fallback (its only use of them).
#
function _positions_from_atoms(atoms::Vector{Atom{Nothing}}, ::Val{D}, ::Type{T}) where {D,T}
    D == 3 && return SVector{3,T}[SVector{3,T}(a.x, a.y, a.z) for a in atoms]
    return SVector{2,T}[SVector{2,T}(a.x, a.y) for a in atoms]
end

"""
    write_output(packmol_system::PackmolSystem)

Write packed coordinates to the output file specified in `packmol_system.output_file`.
Residue numbering follows each structure type's `residue_numbering` (`resnumbers`)
setting, chain identifiers follow `chain`/`changechains`, and a `CRYST1` record is
added when `add_box_sides` is set or the system has periodic boundary conditions.
Also writes any configured `restart_to` file(s) from the current molecule positions,
unless `write_restart=false`.
Returns the path of the file that was written. See [`get_atoms`](@ref) for how the
written coordinates are computed (including PBC wrapping and
`packmol_system.periodic_boundary_style`).
"""
function write_output(packmol_system::PackmolSystem{D,T}; output_file=packmol_system.output_file, write_restart::Bool=true) where {D,T}
    output_atoms = get_atoms(packmol_system)
    write_restart && _write_restart_files(packmol_system, output_atoms)
    atom_positions = _positions_from_atoms(output_atoms, Val(D), T)

    if _file_format(output_file, packmol_system.filetype) == "xyz"
        write_xyz(output_file, Structure(output_atoms; title="Built with Packmol"))
        return output_file
    end
    cryst1_lines = _cryst1_lines(packmol_system, atom_positions)
    header = join(vcat("HEADER", "REMARK   Packmol.jl generated PDB file", cryst1_lines), "\n")
    connectivity = _system_connectivity(packmol_system)
    # Positional `MetaData` (unit cell, connectivity, title): a keyword call with the
    # `Union`-typed connectivity is not compatible with `juliac --trim`
    write_pdb(output_file, Structure(output_atoms, PDBToolsBase.MetaData(nothing, connectivity, nothing)); header)
    if packmol_system.add_amber_ter || packmol_system.amber_ter_preserve
        _add_ter_records(output_file, packmol_system)
    end
    return output_file
end

# Number of copies of a structure type in the output
_ncopies(st::StructureType) = st.fixed.fixed ? 1 : st.number_of_molecules

#
# Connectivity of the packed system (in the order of the atoms of `get_atoms`):
# the connectivity of the template of each structure type, repeated for every
# copy of the molecule, with the indices shifted to the atoms of the copy.
# `nothing` if no structure type has connectivity, or if `connect` is disabled.
#
function _system_connectivity(packmol_system::PackmolSystem)
    packmol_system.connect || return nothing
    any(st -> !isnothing(st.connectivity), packmol_system.structure_types) || return nothing
    connectivity = Vector{Int32}[]
    offset = Int32(0)
    for st in packmol_system.structure_types
        template = st.connectivity
        for _ in 1:_ncopies(st)
            for iat in 1:st.natoms
                bonded = isnothing(template) || iat > length(template) ? Int32[] : template[iat] .+ offset
                push!(connectivity, bonded)
            end
            offset += Int32(st.natoms)
        end
    end
    return connectivity
end

#
# Number of atoms preceding each `TER` record of a structure's PDB file (so
# that 0 is a `TER` before its first atom), for `amber_ter_preserve`.
#
function _ter_positions(filename::AbstractString)
    positions = Int[]
    isfile(filename) || return positions
    natoms = 0
    for line in eachline(filename)
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            natoms += 1
        elseif startswith(line, "TER")
            push!(positions, natoms)
        end
    end
    return positions
end

#
# Add `TER` records, as Fortran Packmol, to the PDB file written by `write_pdb`:
# after every molecule with `add_amber_ter`, and where the structure's own PDB
# file has them with `amber_ter_preserve`.
#
function _add_ter_records(output_file::AbstractString, packmol_system::PackmolSystem)
    # nter[k+1]: number of TER records after the k-th atom of the file
    nter = zeros(Int, sum(st -> st.natoms * _ncopies(st), packmol_system.structure_types; init=0) + 1)
    offset = 0
    for st in packmol_system.structure_types
        ter_positions = packmol_system.amber_ter_preserve ? _ter_positions(st.filename) : Int[]
        for _ in 1:_ncopies(st)
            for k in ter_positions
                k <= st.natoms && (nter[offset+k+1] += 1)
            end
            offset += st.natoms
            packmol_system.add_amber_ter && (nter[offset+1] += 1)
        end
    end
    lines = open(_read_lines, output_file)
    io = open(output_file, "w")
    try
        iat = 0
        for line in lines
            is_atom = startswith(line, "ATOM") || startswith(line, "HETATM")
            if is_atom && iat == 0
                for _ in 1:nter[1]
                    println(io, "TER")
                end
            end
            println(io, line)
            if is_atom
                iat += 1
                for _ in 1:nter[iat+1]
                    println(io, "TER")
                end
            end
        end
    finally
        close(io)
    end
    return output_file
end
