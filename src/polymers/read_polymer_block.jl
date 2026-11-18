#
# The polymers of the input files: the `repeat` blocks of a structure block, and the
# expansion of the polymer into its monomers (one structure type each).
#
#     structure polymer.pdb      # the polymer (one model per conformation)
#       number 10                # polymers
#       repeat 20                # the sequence: repeat blocks of monomers
#         monomer 18             # atoms of the monomer, in the order of the structure
#           link 18 1            # atom 18 of this monomer to atom 1 of the next
#           link 17 2 3.8 10.0   # with the distance and the constant of the link (and,
#                                # optionally, the half-width of a flat bottom)
#       end repeat
#       cyclic                   # the last monomer linked to the first one
#     end structure
#
# The links of a monomer join it to the next monomer of the sequence (of the same repeat
# block or of the next one); those of the last monomer of a linear polymer are ignored.
#

const _default_link_constant = 10.0

#
# Removes the lines of the sequence of a polymer (`repeat`, `monomer`, `link`,
# `end repeat`, `cyclic`) from a structure block, returning the other lines and the
# sequence (empty if the block has no `repeat` blocks).
#
function _read_polymer_spec(input_file_block::IOBuffer, ::Type{T}) where {T}
    seekstart(input_file_block)
    lines = _read_lines(input_file_block)
    filtered = IOBuffer()
    monomers = MonomerSpec{T}[]
    cyclic = false
    in_repeat = false
    nrepeat = 0
    block = MonomerSpec{T}[]
    for line in lines
        words = split(_strip_comment(line))
        if isempty(words)
            println(filtered, line)
            continue
        end
        keyword = words[1]
        if keyword == "repeat"
            in_repeat && throw(ArgumentError("""\n
                repeat blocks cannot be nested: close the previous one with `end repeat`.

            """))
            length(words) >= 2 || throw(ArgumentError("""\n
                repeat requires the number of repeats: $(strip(line))

            """))
            nrepeat = _parse_value(Int, "repeat", words[2]; _val_check=x -> x >= 1 || throw(ArgumentError("""\n
                The number of repeats must be at least 1: $(strip(line))

            """)))
            in_repeat = true
            empty!(block)
        elseif keyword == "end" && length(words) >= 2 && words[2] == "repeat"
            in_repeat || throw(ArgumentError("""\n
                `end repeat` without `repeat`.

            """))
            isempty(block) && throw(ArgumentError("""\n
                A repeat block without monomers.

            """))
            for _ in 1:nrepeat
                for m in block
                    push!(monomers, m)
                end
            end
            in_repeat = false
        elseif keyword == "monomer"
            in_repeat || throw(ArgumentError("""\n
                `monomer` must be inside a repeat block: $(strip(line))

            """))
            length(words) >= 2 || throw(ArgumentError("""\n
                monomer requires its number of atoms: $(strip(line))

            """))
            n = _parse_value(Int, "monomer", words[2]; _val_check=x -> x >= 1 || throw(ArgumentError("""\n
                A monomer must have at least one atom: $(strip(line))

            """)))
            push!(block, MonomerSpec{T}(n, PolymerLink{T}[]))
        elseif keyword == "link"
            (in_repeat && !isempty(block)) || throw(ArgumentError("""\n
                `link` must follow a `monomer` line, inside a repeat block: $(strip(line))

            """))
            (3 <= length(words) <= 6) || throw(ArgumentError("""\n
                link requires two atoms, and optionally the distance, the constant and the
                half-width of the flat bottom: $(strip(line))

            """))
            a = _parse_value(Int, "link", words[2])
            b = _parse_value(Int, "link", words[3])
            d0 = length(words) >= 4 ? _parse_value(T, "link", words[4]) : T(NaN)
            k = length(words) >= 5 ? _parse_value(T, "link", words[5]) : T(NaN)
            w = length(words) >= 6 ? _parse_value(T, "link", words[6]) : T(NaN)
            (1 <= a <= block[end].natoms) || throw(ArgumentError("""\n
                link: atom $a is not an atom of the monomer, of $(block[end].natoms) atoms: $(strip(line))

            """))
            push!(block[end].links, PolymerLink{T}(a, b, d0, k, w))
        elseif keyword == "cyclic"
            cyclic = true
        else
            in_repeat && throw(ArgumentError("""\n
                Keyword $keyword inside a repeat block (only `monomer` and `link`): $(strip(line))

            """))
            println(filtered, line)
        end
    end
    in_repeat && throw(ArgumentError("""\n
        repeat block not closed with `end repeat`.

    """))
    cyclic && isempty(monomers) && throw(ArgumentError("""\n
        cyclic requires the sequence of the polymer (repeat blocks).

    """))
    # Each repeat pushed copies of the same `MonomerSpec`s: independent copies
    monomers = [MonomerSpec{T}(m.natoms, copy(m.links)) for m in monomers]
    return seekstart(filtered), PolymerSpec{T}(monomers, cyclic)
end

#
# The structure types of a structure block of an input file, added to the system: the
# block of a polymer adds its monomers, one structure type each, and the polymer.
#
function _add_structure_block!(sys::PackmolSystem{D,T}, input_file_block::IOBuffer; path::String="") where {D,T}
    block, spec = _read_polymer_spec(input_file_block, T)
    types = _read_structure_types(block, sys.tolerance, Val(D), T; path, filetype=sys.filetype)
    if isempty(spec.monomers)
        append!(sys.structure_types, types)
    else
        _add_polymer!(sys, types, spec)
    end
    return sys
end

#
# Adds a polymer to the system: `rigid_types` are the structure types of the whole
# polymer (one per conformation), and `spec` its sequence of monomers. The monomers are
# new structure types, with the geometry of the first conformation; the distances of the
# links not given are those of the first conformation.
#
function _add_polymer!(sys::PackmolSystem{D,T}, rigid_types::Vector{StructureType{D,T}}, spec::PolymerSpec{T}) where {D,T}
    st = first(rigid_types)
    name = basename(st.filename)
    D == 3 || throw(ArgumentError("""\n
        Polymers ($name) are supported in three dimensions only.

    """))
    st.fixed.fixed && throw(ArgumentError("""\n
        A polymer ($name) cannot be fixed.

    """))
    isnothing(st.restart_from) || throw(ArgumentError("""\n
        restart_from is not supported for polymers ($name).

    """))
    isempty(st.rotation_constraints) || throw(ArgumentError("""\n
        constrain_rotation is not supported for polymers ($name).

    """))
    nmon = length(spec.monomers)
    natoms = sum(m.natoms for m in spec.monomers)
    natoms == st.natoms || throw(ArgumentError("""\n
        The monomers of the polymer $name have $natoms atoms, and its structure has $(st.natoms).

    """))
    nmon >= 2 || throw(ArgumentError("""\n
        The polymer $name must have at least two monomers.

    """))
    # The atoms of each monomer, and the monomer of each atom
    monomer_atoms = UnitRange{Int}[]
    first_atom = 1
    for m in spec.monomers
        push!(monomer_atoms, first_atom:first_atom+m.natoms-1)
        first_atom += m.natoms
    end
    monomer_of_atom = zeros(Int, natoms)
    for (m, r) in enumerate(monomer_atoms)
        monomer_of_atom[r] .= m
    end
    # The links of each junction in each conformation (the distances not given are those
    # of the conformation, not yet centered)
    links = Vector{Vector{PolymerLink{T}}}[]
    for rst in rigid_types
        x = rst.reference_coordinates
        conformation_links = Vector{PolymerLink{T}}[]
        for m in 1:nmon
            next = m < nmon ? m + 1 : (spec.cyclic ? 1 : 0)
            resolved = PolymerLink{T}[]
            if next > 0
                for l in spec.monomers[m].links
                    # (few interpolations: a long string is not compiled by `juliac --trim`)
                    if !(1 <= l.b <= spec.monomers[next].natoms)
                        where = string("link ", l.a, " ", l.b, " of monomer ", m, " of ", name)
                        throw(ArgumentError("""\n
                            $where: atom $(l.b) is not an atom of the next monomer, of $(spec.monomers[next].natoms) atoms.

                        """))
                    end
                    ia = monomer_atoms[m][l.a]
                    ib = monomer_atoms[next][l.b]
                    d0 = isnan(l.d0) ? norm(x[ib] - x[ia]) : l.d0
                    k = isnan(l.k) ? T(_default_link_constant) : l.k
                    w = isnan(l.w) ? zero(T) : l.w
                    push!(resolved, PolymerLink{T}(l.a, l.b, d0, k, w))
                end
            end
            push!(conformation_links, resolved)
        end
        push!(links, conformation_links)
    end
    all(isempty, first(links)) && throw(ArgumentError("""\n
        The polymer $name has no links between its monomers.

    """))
    # The pairs of atoms of different monomers (up to 3 monomers apart) closer than the
    # sum of their radii, in each conformation
    exclusions = Set{Tuple{Int32,Int32}}[]
    excluded_atom = Vector{Bool}[]
    for rst in rigid_types
        x = rst.reference_coordinates
        conformation_exclusions = Set{Tuple{Int32,Int32}}()
        conformation_excluded_atom = zeros(Bool, natoms)
        push!(exclusions, conformation_exclusions)
        push!(excluded_atom, conformation_excluded_atom)
        for i in 1:natoms, j in i+1:natoms
            mi, mj = monomer_of_atom[i], monomer_of_atom[j]
            mi == mj && continue
            sep = abs(mi - mj)
            spec.cyclic && (sep = min(sep, nmon - sep))
            sep <= 3 || continue
            if norm(x[i] - x[j]) < st.radii[i] + st.radii[j]
                push!(conformation_exclusions, (Int32(i), Int32(j)))
                conformation_excluded_atom[i] = conformation_excluded_atom[j] = true
            end
        end
    end
    # The monomers: one structure type for each monomer of each conformation, with a
    # molecule for each polymer of the conformation
    first_type = length(sys.structure_types) + 1
    for rst in rigid_types
        x = rst.reference_coordinates
        for r in monomer_atoms
            atoms = [copy(rst.atoms[i]) for i in r]
            push!(sys.structure_types, StructureType{D,T}(;
                filename=rst.filename, natoms=length(r), atoms,
                number_of_molecules=rst.number_of_molecules,
                reference_coordinates=x[r],
                radii=st.radii[r],
                # the residue numbers of the structure of the polymer
                residue_numbering=st.residue_numbering == -1 ? 1 : st.residue_numbering,
                chain=st.chain, changechains=st.changechains, segid=st.segid,
                # no inflated radii: they would push consecutive monomers apart, against the links
                radscale=something(st.radscale, one(T)),
                connectivity=nothing,
                constraints=copy(st.constraints),
                atom_constraints=[copy(st.atom_constraints[i]) for i in r],
                move_bad_molecules=st.move_bad_molecules,
            ))
        end
    end
    push!(sys.polymers, Polymer{D,T}(
        rigid_types, first_type, monomer_atoms, links, spec.cyclic, exclusions, excluded_atom,
    ))
    return sys
end
