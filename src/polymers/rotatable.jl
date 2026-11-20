#
# Molecules with rotatable bonds (`rotatable i j ...` in a structure block, or
# `flexible_structure` in Julia): the bonds of the molecule (its CONECT records, or the
# pairs of atoms closer than `bond_distance`) are cut at the rotatable bonds, and the
# connected parts of what is left are the rigid fragments of the molecule. Each rotatable
# bond i–j is kept by links: i–j (the bond), and the bonds of i and of j to their other
# neighbors in their fragments, i–m and n–j (the angles of the bond), such that the
# fragments turn about the bond only. The pairs of atoms of different fragments up to
# three bonds apart (1–2, 1–3 and 1–4 pairs) are excluded from the distances.
#

# Whether an atom is a hydrogen: by its element, or else by its name (H, or a digit and H)
function _is_hydrogen(at::Atom)
    element = strip(at.pdb_element)
    isempty(element) || return uppercase(element) == "H"
    name = lstrip(at.name, ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'))
    return startswith(name, "H")
end

# The bonds of a molecule: for each atom, the atoms bonded to it
function _molecule_bonds(st::StructureType{D,T}, bond_distance::T) where {D,T}
    n = st.natoms
    name = basename(st.filename)
    bonded = [Int[] for _ in 1:n]
    if !isnan(bond_distance)
        # the pairs of atoms closer than `bond_distance`, in the first conformation, except
        # pairs of hydrogens (the hydrogens of a CH2 are 1.75 Å apart, as close as the
        # longest bonds of heavy atoms, P–O of 1.6 Å)
        x = st.reference_coordinates
        hydrogen = [_is_hydrogen(at) for at in st.atoms]
        for i in 1:n, j in i+1:n
            (hydrogen[i] && hydrogen[j]) && continue
            if norm(x[i] - x[j]) <= bond_distance
                push!(bonded[i], j)
                push!(bonded[j], i)
            end
        end
    else
        connectivity = st.connectivity
        isnothing(connectivity) && throw(ArgumentError("""\n
            The rotatable bonds of $name require its bonds: the structure file has no CONECT
            records (or `connect no` was set). Give `bond_distance <d>` to compute them from the
            distances (atoms closer than d Å are bonded).

        """))
        for i in 1:n, j in connectivity[i]
            (1 <= j <= n && j != i) || continue
            j in bonded[i] || push!(bonded[i], j)
            i in bonded[j] || push!(bonded[j], i)
        end
    end
    return bonded
end

# The element of an atom (upper case; a String, for `juliac --trim`): by its element, or
# else by the first letter of its name after the digits ("" if none)
function _element(at::Atom)
    element = strip(at.pdb_element)
    isempty(element) || return String(uppercase(element))
    name = lstrip(at.name, ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'))
    return isempty(name) ? "" : string(uppercase(first(name)))
end

#
# The rotatable bonds of `rotatable all`: the bonds not in rings with another heavy atom
# bonded to each end (the bonds that turn more than hydrogens), and certainly single,
# from the numbers of bonds of their atoms: carbon atoms with fewer than four bonds and
# nitrogen atoms with fewer than four bonds are unsaturated (or may be: an amine nitrogen
# has three bonds, as an amide nitrogen), and a bond between two unsaturated atoms does
# not turn (double and triple bonds, amides, and the single bonds that may be double, or
# are conjugated: aryl-vinyl, aryl-carbonyl), unless both are in aromatic rings (rings of
# unsaturated atoms: biaryls). This requires all the hydrogens: without them, a CH2 carbon
# has two bonds (an error if the molecule has no hydrogens).
#
function _all_rotatable_bonds(st::StructureType, bonded::Vector{Vector{Int}})
    n = length(bonded)
    heavy = [!_is_hydrogen(at) for at in st.atoms]
    name = basename(st.filename)
    all(heavy) && throw(ArgumentError("""\n
        rotatable all: $name has no hydrogen atoms. The rotatable bonds are found from the
        numbers of bonds of the atoms, which requires all the hydrogens of the molecule.

    """))
    element = String[_element(at) for at in st.atoms]
    # the bonds in rings (j reachable from i without the bond i-j), and the atoms in rings
    seen = falses(n)
    stack = Int[]
    ring_bond = Set{Tuple{Int,Int}}()
    for i in 1:n, j in bonded[i]
        j > i || continue
        fill!(seen, false)
        seen[i] = true
        empty!(stack)
        push!(stack, i)
        found = false
        while !isempty(stack) && !found
            a = pop!(stack)
            for c in bonded[a]
                (a == i && c == j) && continue
                if c == j
                    found = true
                    break
                end
                seen[c] || (seen[c] = true; push!(stack, c))
            end
        end
        found && push!(ring_bond, (i, j))
    end
    in_ring = falses(n)
    for (i, j) in ring_bond
        in_ring[i] = true
        in_ring[j] = true
    end
    unsaturated(i) = (element[i] == "C" || element[i] == "N") && length(bonded[i]) < 4
    # in an aromatic ring: unsaturated, and bonded in a ring to unsaturated atoms only
    function aromatic(i)
        (unsaturated(i) && in_ring[i]) || return false
        for k in bonded[i]
            ((min(i, k), max(i, k)) in ring_bond && !unsaturated(k)) && return false
        end
        return true
    end
    pairs = Tuple{Int,Int}[]
    for i in 1:n, j in bonded[i]
        j > i || continue
        (heavy[i] && heavy[j]) || continue
        (i, j) in ring_bond && continue
        any(k -> k != j && heavy[k], bonded[i]) || continue
        any(k -> k != i && heavy[k], bonded[j]) || continue
        (unsaturated(i) && unsaturated(j) && !(aromatic(i) && aromatic(j))) && continue
        push!(pairs, (i, j))
    end
    return pairs
end

#
# The fragments of a molecule with its bonds `bonded` cut at the `rotatable` bonds: the
# fragment of each atom, and the atoms of each fragment (in the order of the atoms,
# fragments ordered by their first atom)
#
function _fragments(bonded::Vector{Vector{Int}}, rotatable::Vector{Tuple{Int,Int}})
    n = length(bonded)
    is_cut(i, j) = (i, j) in rotatable || (j, i) in rotatable
    fragment_of = zeros(Int, n)
    fragments = Vector{Int}[]
    for start in 1:n
        fragment_of[start] == 0 || continue
        push!(fragments, Int[])
        f = length(fragments)
        fragment_of[start] = f
        stack = [start]
        while !isempty(stack)
            i = pop!(stack)
            push!(fragments[f], i)
            for j in bonded[i]
                (fragment_of[j] == 0 && !is_cut(i, j)) || continue
                fragment_of[j] = f
                push!(stack, j)
            end
        end
        sort!(fragments[f])
    end
    return fragment_of, fragments
end

#
# The pairs of atoms (i < j) of different fragments up to three bonds apart, given the
# bonds of the molecule and the fragment of each atom
#
function _bonded_exclusions(bonded::Vector{Vector{Int}}, fragment_of::Vector{Int})
    exclusions = Set{Tuple{Int32,Int32}}()
    for i in eachindex(bonded)
        # the atoms within three bonds of i (breadth-first)
        depth = Dict(i => 0)
        queue = [i]
        while !isempty(queue)
            a = popfirst!(queue)
            depth[a] == 3 && continue
            for b in bonded[a]
                haskey(depth, b) && continue
                depth[b] = depth[a] + 1
                push!(queue, b)
            end
        end
        for j in keys(depth)
            (j > i && fragment_of[j] != fragment_of[i]) && push!(exclusions, (Int32(i), Int32(j)))
        end
    end
    return exclusions
end

#
# Why two atoms are not bonded, for the error messages: the source of the bonds, and the
# distance of the atoms (in the first conformation)
#
function _not_bonded_reason(st::StructureType, bond_distance, i::Int, j::Int)
    d = string(round(norm(st.reference_coordinates[i] - st.reference_coordinates[j]); digits=2))
    isnan(bond_distance) && return string("not bonded in the CONECT records (they are ", d, " Å apart)")
    return string("farther than bond_distance ", bond_distance, " Å (they are ", d, " Å apart)")
end

#
# The pairs of atoms of the links of a rotatable bond i-j between two fragments: the bond,
# and the bonds of each end to its other neighbors in its fragment (which keep the angles
# of the bond)
#
function _rotatable_link_pairs!(pairs::Vector{Tuple{Int,Int}}, i::Int, j::Int, bonded::Vector{Vector{Int}}, fragment_of::Vector{Int})
    push!(pairs, (i, j))
    for m in bonded[i]
        (m != j && fragment_of[m] == fragment_of[i]) && push!(pairs, (m, j))
    end
    for m in bonded[j]
        (m != i && fragment_of[m] == fragment_of[j]) && push!(pairs, (i, m))
    end
    return pairs
end

#
# Adds a molecule with rotatable bonds to the system: `rigid_types` are the structure
# types of its conformations
#
function _add_rotatable!(
    sys::PackmolSystem{D,T}, rigid_types::Vector{StructureType{D,T}}, rotatable::Vector{Tuple{Int,Int}}, bond_distance::T,
) where {D,T}
    st = first(rigid_types)
    name = basename(st.filename)
    _check_flexible(st)
    n = st.natoms
    bonded = _molecule_bonds(st, bond_distance)
    if (0, 0) in rotatable
        # `rotatable all`: the bonds found, and the others given
        found = _all_rotatable_bonds(st, bonded)
        isempty(found) && throw(ArgumentError("""\n
            rotatable all: no rotatable bonds found in $name.

        """))
        given = filter(b -> b != (0, 0) && !((b[1], b[2]) in found || (b[2], b[1]) in found), rotatable)
        rotatable = vcat(found, given)
        _println("  Rotatable bonds of ", name, ": ", length(found), " found")
    end
    for (i, j) in rotatable
        (1 <= i <= n && 1 <= j <= n) || throw(ArgumentError("""\n
            $(string("rotatable ", i, " ", j)): the molecule $name has $n atoms.

        """))
        # (few interpolations: a long string is not compiled by `juliac --trim`)
        bond = string("rotatable ", i, " ", j)
        reason = j in bonded[i] ? "" : _not_bonded_reason(st, bond_distance, i, j)
        isempty(reason) || throw(ArgumentError("""\n
            $bond: the atoms are $reason, in $name.

        """))
    end
    fragment_of, fragments = _fragments(bonded, rotatable)
    for (i, j) in rotatable
        bond = string("rotatable ", i, " ", j)
        fragment_of[i] == fragment_of[j] && throw(ArgumentError("""\n
            $bond: the bond is in a ring of $name (its atoms are in the same rigid fragment
            once the rotatable bonds are cut), and cannot rotate.

        """))
    end
    # The index of each atom in its fragment
    local_index = zeros(Int, n)
    for r in fragments, (k, i) in enumerate(r)
        local_index[i] = k
    end
    # The pairs of atoms of each link: the bond, and the bonds of each end to its other
    # neighbors in its fragment
    pairs = Tuple{Int,Int}[]
    for (i, j) in rotatable
        _rotatable_link_pairs!(pairs, i, j, bonded, fragment_of)
    end
    links = Vector{FragmentLink{T}}[]
    for rst in rigid_types
        x = rst.reference_coordinates
        push!(links, [FragmentLink{T}(fragment_of[a], local_index[a], fragment_of[b], local_index[b],
            norm(x[b] - x[a]), T(_default_link_constant), zero(T)) for (a, b) in pairs])
    end
    exclusions = _bonded_exclusions(bonded, fragment_of)
    return _add_flexible_molecule!(sys, rigid_types, fragments, links, fill(exclusions, length(rigid_types)), false)
end
