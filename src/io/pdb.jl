#
# PDB reading and writing, adapted from PDBTools.jl (`read_pdb`/`write_pdb`) to what
# Packmol needs: plain `ATOM`/`HETATM` records into `Atom{Nothing}` vectors, without
# selections or mmCIF support. Kept type-stable so that it compiles with `juliac --trim`.
#
# Differences from PDBTools:
#  - two consecutive atoms are in the same residue if their residue numbers, residue
#    names, chains, segment names and models are the same (PDBTools also merges the
#    alternate conformations of protein residues);
#  - the element column is written as read from the input file (blank if it was not
#    there), instead of being inferred from the atom name.
#

const _PDBString = fieldtype(Atom{Nothing}, :name)

# The `range` columns of `line`, or the part of them that the line contains.
function _columns(line::AbstractString, range::UnitRange{Int})
    first(range) > lastindex(line) && return SubString(line, 1, 0)
    return SubString(line, first(range), min(last(range), lastindex(line)))
end

function _pdb_string(line::AbstractString, range::UnitRange{Int})
    s = strip(_columns(line, range))
    return isempty(s) ? _PDBString("X") : _PDBString(s)
end

# Integers might be written in hexadecimal when they don't fit the columns.
function _pdb_int(line::AbstractString, range::UnitRange{Int})
    s = strip(_columns(line, range))
    i = tryparse(Int32, s)
    isnothing(i) && (i = tryparse(Int32, s; base=16))
    return something(i, Int32(0))
end

function _pdb_float(line::AbstractString, range::UnitRange{Int}; required::Bool=false)
    s = strip(_columns(line, range))
    x = tryparse(Float32, s)
    if isnothing(x)
        required && throw(ArgumentError("Could not read coordinate from PDB line: \"$line\""))
        return 0.0f0
    end
    return x
end

# The charge sometimes appears as "+1" or "-1" and sometimes as "1+" or "1-".
function _pdb_charge(line::AbstractString, range::UnitRange{Int})
    s = strip(_columns(line, range))
    sign = occursin('-', s) ? -1.0f0 : 1.0f0
    x = tryparse(Float32, filter(c -> c != '+' && c != '-', s))
    return isnothing(x) ? 0.0f0 : sign * x
end

# The position of an atom, as an `SVector{D,T}`
_atom_position(atom::Atom, ::Val{D}, ::Type{T}) where {D,T} = SVector{D,T}(ntuple(i -> (atom.x, atom.y, atom.z)[i], D))

function _same_residue(at1::Atom{Nothing}, at2::Atom{Nothing})
    at1.resnum == at2.resnum && at1.resname == at2.resname && at1.chain == at2.chain &&
        at1.segname == at2.segname && at1.model == at2.model
end

"""
    read_pdb_atoms(filename::AbstractString)

Read the `ATOM` and `HETATM` records of the PDB file `filename` into a `Vector{Atom{Nothing}}`.

"""
function read_pdb_atoms(filename::AbstractString)
    atoms = Atom{Nothing}[]
    model = Int32(1)
    for line in eachline(filename)
        if occursin("END", line)
            model += Int32(1)
        end
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            atom = Atom{Nothing}(
                Int32(length(atoms) + 1),         # index
                _pdb_int(line, 7:11),             # index_pdb
                _pdb_string(line, 13:16),         # name
                _pdb_string(line, 17:21),         # resname
                _pdb_string(line, 22:22),         # chain
                _pdb_int(line, 23:26),            # resnum
                Int32(1),                         # residue
                _pdb_float(line, 31:38; required=true), # x
                _pdb_float(line, 39:46; required=true), # y
                _pdb_float(line, 47:54; required=true), # z
                _pdb_float(line, 61:66),          # beta
                _pdb_float(line, 55:60),          # occup
                model,                            # model
                _pdb_string(line, 73:76),         # segname
                _pdb_string(line, 77:78),         # pdb_element
                _pdb_charge(line, 79:80),         # charge
                nothing,                          # custom
                Int8(startswith(line, "HETATM")), # flag
            )
            if !isempty(atoms)
                last_atom = atoms[end]
                atom.residue = _same_residue(atom, last_atom) ? last_atom.residue : last_atom.residue + Int32(1)
            end
            push!(atoms, atom)
        end
    end
    if isempty(atoms)
        throw(ArgumentError("""\n
            Could not find any atom in PDB file $filename.

        """))
    end
    return atoms
end

function _align_name(name::AbstractString)
    name = strip(name)
    length(name) == 1 && return " $(name)  "
    length(name) == 2 && return " $(name) "
    length(name) == 3 && return " $(name)"
    return String(name)
end

function _align_resname(resname::AbstractString)
    resname = strip(resname)
    length(resname) == 1 && return " $(resname) "
    length(resname) == 2 && return "$(resname)  "
    length(resname) == 3 && return "$(resname) "
    return String(resname)
end

#
# The `ATOM`/`HETATM` line of an atom in PDB format. Indices and residue numbers that
# don't fit the columns are written in hexadecimal.
#
#COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#73 - 76        String        segname      Segment identifier, left-justified (not default of PDB)
#77 - 78        LString(2)    element      Element symbol string, right-justified.
#
function pdb_atom_line(atom::Atom{Nothing})
    record = atom.flag == 0 ? "ATOM" : "HETATM"
    index = mod(atom.index, 1048576)
    name = _align_name(atom.name)
    resname = _align_resname(atom.resname)
    chain = String(atom.chain)
    resnum = mod(atom.resnum, 65535)
    segname = String(atom.segname)
    element = atom.pdb_element == "X" ? "" : String(atom.pdb_element)
    index_str = atom.index <= 99999 ? @sprintf("%5i", index) : @sprintf("%5x", index)
    resnum_str = atom.resnum <= 9999 ? @sprintf("%4i", resnum) : @sprintf("%4x", resnum)
    return @sprintf(
        "%-6s%5s%1s%-4s%1s%4s%1s%4s%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s",
        record, index_str, " ", name, " ", resname, chain, resnum_str, "    ",
        atom.x, atom.y, atom.z, atom.occup, atom.beta, "      ", segname, element,
    )
end

"""
    write_pdb_atoms(filename::AbstractString, atoms::AbstractVector{Atom{Nothing}}; header::String="")

Write `atoms` to the PDB file `filename`, preceded by the `header` lines, if any, and
followed by an `END` line.

"""
function write_pdb_atoms(filename::AbstractString, atoms::AbstractVector{Atom{Nothing}}; header::String="")
    open(filename; write=true) do io
        isempty(header) || println(io, header)
        for atom in atoms
            println(io, pdb_atom_line(atom))
        end
        println(io, "END")
    end
    return filename
end
