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
# Restart file I/O (Fortran Packmol's `restart_from`/`restart_to`): a plain
# text dump, one line per molecule, of that molecule's raw optimization
# variables (center of mass + Euler angles, in radians) — not a PDB file.
# Written every time `write_output` runs (matching Fortran, which writes
# restart files at the top of its own `output` subroutine, called for every
# intermediate/best-so-far write as well as the final one), and read back
# in `packmol()` before the initial-approximation pipeline runs, so a
# restart resumes packing from exactly where the previous run left off
# instead of re-deriving a starting point.
#
function _write_restart(filename::AbstractString, positions::AbstractVector{MoleculePosition{D,T}}) where {D,T}
    open(filename, "w") do io
        for mp in positions
            for v in mp.cm
                @printf(io, " %23.16e", v)
            end
            for v in mp.angles
                @printf(io, " %23.16e", v)
            end
            println(io)
        end
    end
    return filename
end

function _read_restart(filename::AbstractString, n::Int, ::Type{MoleculePosition{D,T}}) where {D,T}
    positions = Vector{MoleculePosition{D,T}}(undef, n)
    open(filename) do io
        for i in 1:n
            eof(io) && error("Restart file '$filename' has fewer than the expected $n molecule lines.")
            vals = parse.(T, split(readline(io)))
            length(vals) < 2D && error("Restart file '$filename' line $i has fewer than $(2D) values.")
            positions[i] = MoleculePosition(SVector{D,T}(vals[1:D]), SVector{D,T}(vals[D+1:2D]))
        end
    end
    return positions
end

#
# Recovers a molecule's (cm, angles) from its atoms' observed Cartesian
# positions, by rigid-body-aligning `reference` (that structure type's own
# `reference_coordinates` — already centered at its own centroid) onto
# `observed` (raw, uncentered) via the Kabsch algorithm: `cm` is `observed`'s
# centroid, and `angles` is the least-squares-optimal rotation taking
# `reference` onto `observed - cm`, decomposed back into `eulermat`'s
# parameterization. Used for a PDB-based restart_from, where the exact
# rigid-body transform used to generate the observed positions is unknown
# and must be recovered rather than read out directly (unlike the raw
# restart format, which stores it verbatim).
#
function _align_molecule(reference::AbstractVector{SVector{3,T}}, observed::AbstractVector{SVector{3,T}}) where {T}
    cm = sum(observed) / length(observed)
    R = kabsch_rotation(reference, observed .- Ref(cm))
    return MoleculePosition(cm, euler_angles(R))
end

#
# Reads restart positions for a sequence of structure-type "segments" — each
# a (number_of_molecules, natoms, reference_coordinates) triple, in the same
# order the atoms appear in `atoms` — by rigid-body-aligning every molecule's
# atoms against its own type's template (see `_align_molecule`). Used for
# both a `.pdb`-file restart_from and an in-memory `Vector{<:Atom}` one; both
# go through this same function once the atoms are in hand.
#
function _restart_positions_from_atoms(
    atoms::AbstractVector{<:Atom}, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}}}},
) where {D,T}
    expected = sum(nmols * natoms for (nmols, natoms, _) in segments)
    length(atoms) == expected || error(
        "Restart PDB has $(length(atoms)) atoms, expected $expected " *
        "($(join(("$nmols × $natoms" for (nmols, natoms, _) in segments), " + ")))."
    )
    positions = MoleculePosition{D,T}[]
    iat = 0
    for (nmols, natoms, refcoords) in segments
        for _ in 1:nmols
            observed = SVector{D,T}[SVector{D,T}(coor(atoms[iat+j])[1:D]) for j in 1:natoms]
            iat += natoms
            push!(positions, _align_molecule(refcoords, observed))
        end
    end
    return positions
end

# A `.pdb` extension (case-insensitive) selects the rigid-body-alignment
# restart path; anything else is read as the raw Packmol restart format.
_is_pdb_path(path::AbstractString) = lowercase(splitext(path)[2]) == ".pdb"

#
# Restart positions for one or more structure-type `segments` (see
# `_restart_positions_from_atoms`), dispatching on what `source` is: an
# in-memory atom vector (Julia API only), a `.pdb` file path (rigid-body
# alignment against each segment's template), or any other file path (the
# raw Packmol restart format, read directly with no alignment needed since
# it already stores (cm, angles) verbatim).
#
function _restart_positions(source::Vector{<:Atom}, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}}}}, ::Type{MoleculePosition{D,T}}) where {D,T}
    _restart_positions_from_atoms(source, segments)
end
function _restart_positions(source::AbstractString, segments::Vector{Tuple{Int,Int,Vector{SVector{D,T}}}}, ::Type{MoleculePosition{D,T}}) where {D,T}
    if _is_pdb_path(source)
        _restart_positions_from_atoms(read_pdb(source), segments)
    else
        _read_restart(source, sum(nmols for (nmols, _, _) in segments), MoleculePosition{D,T})
    end
end

#
# Write any configured restart_to files (whole-system and/or per structure
# type) from the current molecule positions. A no-op for any that aren't set.
#
function _write_restart_files(packmol_system::PackmolSystem{D,T}) where {D,T}
    isnothing(packmol_system.restart_to) ||
        _write_restart(packmol_system.restart_to, packmol_system.molecule_positions)
    imol_offset = 0
    for st in packmol_system.structure_types
        if !isnothing(st.restart_to)
            _write_restart(st.restart_to, view(packmol_system.molecule_positions, imol_offset+1:imol_offset+st.number_of_molecules))
        end
        imol_offset += st.number_of_molecules
    end
    return nothing
end

"""
    write_output(packmol_system::PackmolSystem)

Write packed coordinates to the output file specified in `packmol_system.output_file`.
Residue numbering follows each structure type's `residue_numbering` (`resnumbers`)
setting, chain identifiers follow `chain`/`changechains`, and a `CRYST1` record is
added when `add_box_sides` is set or the system has periodic boundary conditions.
Also writes any configured `restart_to` file(s) from the current molecule positions.
Returns the path of the file that was written.
"""
function write_output(packmol_system::PackmolSystem{D,T}; output_file=packmol_system.output_file) where {D,T}
    _write_restart_files(packmol_system)
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    output_atoms = Atom{Nothing}[]
    iat = 0
    atom_serial = 0
    iimol = 0       # molecule index across the whole system (fixed + free) — resnumbers mode 3
    irescount = 1   # running residue-number offset across the whole system — resnumbers mode 2
    ichain = 0      # count of auto-assigned chain letter groups so far, across all structure types
    for st in packmol_system.structure_types
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
                push!(output_atoms, atom)
            end
            irescount += nres
        else
            odd_chain = even_chain = ' '
            for imol_local in 1:st.number_of_molecules
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
                    push!(output_atoms, atom)
                end
                irescount += nres
            end
        end
    end

    cryst1_lines = _cryst1_lines(packmol_system, atom_positions)
    header = join(vcat("HEADER", "REMARK   Packmol.jl generated PDB file", cryst1_lines), "\n")
    write_pdb(output_file, output_atoms; header)
    return output_file
end
