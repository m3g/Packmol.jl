#
# The flexible molecules of the input files: the `fragment` lines of a structure block
# (alone or in `repeat` blocks), and the expansion of the molecule into its fragments
# (one structure type each). The fragments of a polymer are its monomers.
#
#     structure polymer.pdb      # the polymer (one model per conformation)
#       number 10                # polymers
#       repeat 20                # the sequence: fragments, alone or in repeat blocks
#         fragment 18            # atoms of the fragment, in the order of the structure
#           link 18 1            # atom 18 of this fragment to atom 1 of the next
#           link 17 2 3.8 10.0   # with the distance and the constant of the link (and,
#                                # optionally, the half-width of a flat bottom)
#           rotatable 18 1       # or: the bond of atom 18 to atom 1 of the next fragment,
#                                # with the links that keep its length and angles
#       end repeat
#     (`rotatable i j ...` or `rotatable all`, not after a `fragment` line: the rotatable
#     bonds of a molecule, see rotatable.jl)
#       bond_distance 1.9        # the bonds (for `rotatable`), if no CONECT records
#       cyclic                   # the last fragment linked to the first one
#     end structure
#
# The links of a fragment join it to the next fragment of the sequence (of the same repeat
# block or not); those of the last fragment of a linear molecule are ignored.
#

const _default_link_constant = 10.0

#
# Removes the lines of the fragments of a flexible molecule (`fragment`, `link`,
# `rotatable`, `repeat`, `end repeat`, `cyclic`, `bond_distance`) from a structure block,
# returning the other lines and the specification (empty if the molecule is rigid).
# `fragment` lines, with their `link` and `rotatable` lines, may be inside or outside
# repeat blocks; `rotatable` lines not following a `fragment` line are the rotatable
# bonds of the molecule.
#
function _read_polymer_spec(input_file_block::IOBuffer, ::Type{T}) where {T}
    seekstart(input_file_block)
    lines = _read_lines(input_file_block)
    filtered = IOBuffer()
    fragments = FragmentSpec{T}[]
    cyclic = false
    rotatable = Tuple{Int,Int}[]
    bond_distance = T(NaN)
    in_repeat = false
    nrepeat = 0
    block = FragmentSpec{T}[]
    # the fragment of the last `fragment` line (its `link` and `rotatable` lines follow)
    current = nothing
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
            current = nothing
            empty!(block)
        elseif keyword == "end" && length(words) >= 2 && words[2] == "repeat"
            in_repeat || throw(ArgumentError("""\n
                `end repeat` without `repeat`.

            """))
            isempty(block) && throw(ArgumentError("""\n
                A repeat block without fragments.

            """))
            for _ in 1:nrepeat
                for m in block
                    push!(fragments, m)
                end
            end
            in_repeat = false
            current = nothing
        elseif keyword == "fragment"
            length(words) == 2 || throw(ArgumentError("""\n
                fragment requires its number of atoms: $(strip(line))

            """))
            n = _parse_value(Int, "fragment", words[2]; _val_check=x -> x >= 1 || throw(ArgumentError("""\n
                A fragment must have at least one atom: $(strip(line))

            """)))
            current = FragmentSpec{T}(n, PolymerLink{T}[], Tuple{Int,Int}[])
            push!(in_repeat ? block : fragments, current)
        elseif keyword == "link"
            isnothing(current) && throw(ArgumentError("""\n
                `link` must follow a `fragment` line: $(strip(line))

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
            (1 <= a <= current.natoms) || throw(ArgumentError("""\n
                link: atom $a is not an atom of the fragment, of $(current.natoms) atoms: $(strip(line))

            """))
            push!(current.links, PolymerLink{T}(a, b, d0, k, w))
        elseif keyword == "cyclic"
            cyclic = true
        elseif keyword == "rotatable" && !isnothing(current)
            # the bond of an atom of the fragment to an atom of the next one
            length(words) == 3 || throw(ArgumentError("""\n
                After a `fragment` line, rotatable gives an atom of the fragment and an atom of
                the next one: $(strip(line))

            """))
            a = _parse_value(Int, "rotatable", words[2])
            b = _parse_value(Int, "rotatable", words[3])
            (1 <= a <= current.natoms) || throw(ArgumentError("""\n
                rotatable: atom $a is not an atom of the fragment, of $(current.natoms) atoms: $(strip(line))

            """))
            push!(current.rotatable, (a, b))
        elseif keyword == "rotatable" && in_repeat
            throw(ArgumentError("""\n
                In a repeat block, rotatable follows a `fragment` line: $(strip(line))

            """))
        elseif keyword == "rotatable" && length(words) == 2 && words[2] == "all"
            # every rotatable bond, found from the bonds (`_all_rotatable_bonds`)
            push!(rotatable, (0, 0))
        elseif keyword == "rotatable"
            (length(words) >= 3 && isodd(length(words))) || throw(ArgumentError("""\n
                rotatable requires pairs of bonded atoms: $(strip(line))

            """))
            for k in 2:2:length(words)
                push!(rotatable, (_parse_value(Int, "rotatable", words[k]), _parse_value(Int, "rotatable", words[k+1])))
            end
        elseif keyword == "bond_distance"
            length(words) == 2 || throw(ArgumentError("""\n
                bond_distance requires the largest distance of bonded atoms: $(strip(line))

            """))
            bond_distance = _parse_value(T, "bond_distance", words[2])
        else
            in_repeat && throw(ArgumentError("""\n
                Keyword $keyword inside a repeat block (only `fragment`, `link` and `rotatable`): $(strip(line))

            """))
            current = nothing
            println(filtered, line)
        end
    end
    in_repeat && throw(ArgumentError("""\n
        repeat block not closed with `end repeat`.

    """))
    cyclic && isempty(fragments) && throw(ArgumentError("""\n
        cyclic requires the fragments of the molecule (`fragment` lines).

    """))
    !isempty(rotatable) && !isempty(fragments) && throw(ArgumentError("""\n
        The rotatable bonds of a molecule given by its fragments are given in the fragments
        (`rotatable <a> <b>` after a `fragment` line: atom a of the fragment, atom b of the
        next one).

    """))
    # Each repeat pushed copies of the same `FragmentSpec`s: independent copies
    fragments = [FragmentSpec{T}(m.natoms, copy(m.links), copy(m.rotatable)) for m in fragments]
    return seekstart(filtered), PolymerSpec{T}(fragments, cyclic, rotatable, bond_distance)
end

#
# The structure types of a structure block of an input file, added to the system: the
# block of a flexible molecule (a polymer, or a molecule with rotatable bonds) adds its
# fragments, one structure type each for each conformation, and the molecule.
#
function _add_structure_block!(sys::PackmolSystem{D,T}, input_file_block::IOBuffer; path::String="") where {D,T}
    block, spec = _read_polymer_spec(input_file_block, T)
    types = _read_structure_types(block, sys.tolerance, Val(D), T; path, filetype=sys.filetype)
    _add_flexible!(sys, types, spec)
    return sys
end

# Adds the structure types of a molecule: those of its fragments if it is flexible (a
# polymer, or with rotatable bonds), otherwise its structure types
function _add_flexible!(sys::PackmolSystem{D,T}, types::Vector{StructureType{D,T}}, spec::PolymerSpec{T}) where {D,T}
    if !isempty(spec.fragments)
        _add_polymer!(sys, types, spec)
    elseif !isempty(spec.rotatable)
        _add_rotatable!(sys, types, spec.rotatable, spec.bond_distance)
    else
        append!(sys.structure_types, types)
    end
    return sys
end

# The options not supported for flexible molecules (polymers and molecules with rotatable bonds)
function _check_flexible(st::StructureType{D}) where {D}
    name = basename(st.filename)
    D == 3 || throw(ArgumentError("""\n
        Flexible molecules ($name: polymers, rotatable bonds) are supported in three dimensions only.

    """))
    st.fixed.fixed && throw(ArgumentError("""\n
        A flexible molecule ($name: a polymer, or with rotatable bonds) cannot be fixed.

    """))
    isnothing(st.restart_from) || throw(ArgumentError("""\n
        restart_from is not supported for flexible molecules ($name).

    """))
    isempty(st.rotation_constraints) || throw(ArgumentError("""\n
        constrain_rotation is not supported for flexible molecules ($name).

    """))
    return nothing
end

#
# Adds a polymer to the system: `rigid_types` are the structure types of the whole
# polymer (one per conformation), and `spec` its sequence of monomers. The monomers are
# the fragments of the molecule; the links not given distances have those of each
# conformation; the rotatable bonds of the monomers add the links of each bond (the bond,
# and the bonds of its atoms to their neighbors in their monomers, from the bonds of the
# polymer). The excluded pairs are the atoms of different monomers (up to 3 apart) closer
# than the sum of their radii in a conformation, and, if the bonds are known (rotatable
# bonds), those up to three bonds apart.
#
function _add_polymer!(sys::PackmolSystem{D,T}, rigid_types::Vector{StructureType{D,T}}, spec::PolymerSpec{T}) where {D,T}
    st = first(rigid_types)
    name = basename(st.filename)
    _check_flexible(st)
    nmon = length(spec.fragments)
    natoms = sum(m.natoms for m in spec.fragments)
    natoms == st.natoms || throw(ArgumentError("""\n
        The fragments of $name have $natoms atoms, and its structure has $(st.natoms).

    """))
    nmon >= 2 || throw(ArgumentError("""\n
        The flexible molecule $name must have at least two fragments.

    """))
    # The atoms of each monomer, and the monomer of each atom
    fragments = Vector{Int}[]
    first_atom = 1
    for m in spec.fragments
        push!(fragments, collect(first_atom:first_atom+m.natoms-1))
        first_atom += m.natoms
    end
    monomer_of_atom = zeros(Int, natoms)
    for (m, r) in enumerate(fragments)
        monomer_of_atom[r] .= m
    end
    # The rotatable bonds of the monomers: the pairs of atoms (in the polymer) of their
    # links, from the bonds of the polymer
    has_rotatable = any(m -> !isempty(m.rotatable), spec.fragments)
    bonded = has_rotatable ? _molecule_bonds(st, spec.bond_distance) : Vector{Int}[]
    rotatable_pairs = Tuple{Int,Int}[]
    for m in 1:nmon
        next = m < nmon ? m + 1 : (spec.cyclic ? 1 : 0)
        next > 0 || continue
        for (a, b) in spec.fragments[m].rotatable
            if !(1 <= b <= spec.fragments[next].natoms)
                where = string("rotatable ", a, " ", b, " of fragment ", m, " of ", name)
                throw(ArgumentError("""\n
                    $where: atom $b is not an atom of the next fragment, of $(spec.fragments[next].natoms) atoms.

                """))
            end
            ia = fragments[m][a]
            ib = fragments[next][b]
            if !(ib in bonded[ia])
                where = string("rotatable ", a, " ", b, " of fragment ", m, " of ", name)
                reason = _not_bonded_reason(st, spec.bond_distance, ia, ib)
                throw(ArgumentError("""\n
                    $where (atoms $ia and $ib of the polymer): the atoms are $reason.

                """))
            end
            _rotatable_link_pairs!(rotatable_pairs, ia, ib, bonded, monomer_of_atom)
        end
    end
    # The links of each junction in each conformation (the distances not given are those
    # of the conformation, not yet centered)
    local_index = zeros(Int, natoms)
    for r in fragments, (k, i) in enumerate(r)
        local_index[i] = k
    end
    links = Vector{FragmentLink{T}}[]
    for rst in rigid_types
        x = rst.reference_coordinates
        conformation_links = FragmentLink{T}[]
        for (ia, ib) in rotatable_pairs
            push!(conformation_links, FragmentLink{T}(monomer_of_atom[ia], local_index[ia], monomer_of_atom[ib], local_index[ib],
                norm(x[ib] - x[ia]), T(_default_link_constant), zero(T)))
        end
        for m in 1:nmon
            next = m < nmon ? m + 1 : (spec.cyclic ? 1 : 0)
            next > 0 || continue
            for l in spec.fragments[m].links
                # (few interpolations: a long string is not compiled by `juliac --trim`)
                if !(1 <= l.b <= spec.fragments[next].natoms)
                    where = string("link ", l.a, " ", l.b, " of fragment ", m, " of ", name)
                    throw(ArgumentError("""\n
                        $where: atom $(l.b) is not an atom of the next fragment, of $(spec.fragments[next].natoms) atoms.

                    """))
                end
                ia = fragments[m][l.a]
                ib = fragments[next][l.b]
                d0 = isnan(l.d0) ? norm(x[ib] - x[ia]) : l.d0
                k = isnan(l.k) ? T(_default_link_constant) : l.k
                w = isnan(l.w) ? zero(T) : l.w
                push!(conformation_links, FragmentLink{T}(m, l.a, next, l.b, d0, k, w))
            end
        end
        push!(links, conformation_links)
    end
    isempty(first(links)) && throw(ArgumentError("""\n
        The flexible molecule $name has no links between its fragments.

    """))
    # The pairs of atoms of different monomers (up to 3 monomers apart) closer than the
    # sum of their radii, in each conformation
    exclusions = Set{Tuple{Int32,Int32}}[]
    bonded_exclusions = has_rotatable ? _bonded_exclusions(bonded, monomer_of_atom) : Set{Tuple{Int32,Int32}}()
    for rst in rigid_types
        x = rst.reference_coordinates
        conformation_exclusions = Set{Tuple{Int32,Int32}}()
        for i in 1:natoms, j in i+1:natoms
            mi, mj = monomer_of_atom[i], monomer_of_atom[j]
            mi == mj && continue
            sep = abs(mi - mj)
            spec.cyclic && (sep = min(sep, nmon - sep))
            sep <= 3 || continue
            if norm(x[i] - x[j]) < st.radii[i] + st.radii[j]
                push!(conformation_exclusions, (Int32(i), Int32(j)))
            end
        end
        union!(conformation_exclusions, bonded_exclusions)
        push!(exclusions, conformation_exclusions)
    end
    return _add_flexible_molecule!(sys, rigid_types, fragments, links, exclusions, spec.cyclic)
end

#
# Adds a flexible molecule to the system, from the structure types of its conformations
# (`rigid_types`), its fragments (atoms of the structure), the links between them and
# the excluded pairs (`i < j`) of each conformation: a structure type for each fragment
# of each conformation, with a molecule for each copy of the molecule in that conformation.
#
function _add_flexible_molecule!(
    sys::PackmolSystem{D,T}, rigid_types::Vector{StructureType{D,T}}, fragments::Vector{Vector{Int}},
    links::Vector{Vector{FragmentLink{T}}}, exclusions::Vector{Set{Tuple{Int32,Int32}}}, cyclic::Bool,
) where {D,T}
    st = first(rigid_types)
    excluded_atom = Vector{Bool}[]
    for conformation_exclusions in exclusions
        e = zeros(Bool, st.natoms)
        for (i, j) in conformation_exclusions
            e[i] = e[j] = true
        end
        push!(excluded_atom, e)
    end
    first_type = length(sys.structure_types) + 1
    for rst in rigid_types
        x = rst.reference_coordinates
        for r in fragments
            atoms = [copy(rst.atoms[i]) for i in r]
            push!(sys.structure_types, StructureType{D,T}(;
                filename=rst.filename, natoms=length(r), atoms,
                number_of_molecules=rst.number_of_molecules,
                reference_coordinates=x[r],
                radii=st.radii[r],
                # the residue numbers of the structure of the molecule
                residue_numbering=st.residue_numbering == -1 ? 1 : st.residue_numbering,
                chain=st.chain, changechains=st.changechains, segid=st.segid,
                # no inflated radii: they would push the fragments apart, against the links
                radscale=something(st.radscale, one(T)),
                connectivity=nothing,
                constraints=copy(st.constraints),
                atom_constraints=[copy(st.atom_constraints[i]) for i in r],
                move_bad_molecules=st.move_bad_molecules,
            ))
        end
    end
    push!(sys.polymers, Polymer{D,T}(rigid_types, first_type, fragments, links, cyclic, exclusions, excluded_atom))
    return sys
end
