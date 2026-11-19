#
# The function and gradient of the links between the monomers of the polymers:
# `k (d - d0)²` for each link, added to the function of both molecules. A link deviating
# from its distance by `overlap_precision` (the precision of the tolerance) or more does
# not satisfy the convergence criterion.
#

# The vector from `xa` to `xb` (its shortest periodic image, with PBC)
@inline _link_vector(xa, xb, has_pbc::Bool, unitcell) = has_pbc ? _min_image(xb - xa, unitcell) : xb - xa

# The penalty of a link and the gradient relative to `xb` (that relative to `xa` is its
# opposite), and the deviation of its distance from the flat bottom `d0 ± w`
@inline function _link_penalty(v::SVector{D,T}, l::SystemLink{T}) where {D,T}
    d = norm(v)
    dev = d - l.d0
    excess = abs(dev) - l.w
    excess <= zero(T) && return zero(T), zero(v), zero(T)
    p = l.k * excess^2
    g = d > zero(T) ? (2 * l.k * excess * sign(dev) / d) * v : zero(v)
    return p, g, excess
end

#
# The links of the whole system, with the positions of all the atoms
#
function _link_fg!(fg::InteratomicDistanceFG{D,T}, atom_positions::Vector{SVector{D,T}}, sys::PackmolSystem{D,T}) where {D,T}
    links = sys.polymer_topology.links
    fg.max_link_deviation = zero(T)
    isempty(links) && return fg
    has_pbc, unitcell, _ = _typed_unitcell(sys)
    maxdev = zero(T)
    for l in links
        v = _link_vector(atom_positions[l.ia], atom_positions[l.ib], has_pbc, unitcell)
        p, g, dev = _link_penalty(v, l)
        fg.f += p
        fg.fmol[l.ma] += p
        fg.fmol[l.mb] += p
        fg.gxcar[l.ia] -= g
        fg.gxcar[l.ib] += g
        maxdev = max(maxdev, dev)
        if fg.count_overlaps[] && dev >= fg.overlap_precision
            push!(fg.unconverged_mols, l.ma, l.mb)
        end
    end
    fg.max_link_deviation = maxdev
    return fg
end

#
# A link of the active subsystem (see packing/active_set.jl): the local indices of its
# atoms among the active atoms (0 if frozen), and the positions of the frozen ones
#
struct ActiveLink{D,T}
    link::SystemLink{T}
    ia::Int
    ib::Int
    xa::SVector{D,T}
    xb::SVector{D,T}
end

# The links with an atom of an active molecule
function _active_links(sys::PackmolSystem{D,T}, atom_map::Vector{Int}, cl_positions::Vector{SVector{D,T}}) where {D,T}
    links = sys.polymer_topology.links
    active_links = ActiveLink{D,T}[]
    isempty(links) && return active_links
    local_index = zeros(Int, length(sys.atoms))
    for (il, iat) in enumerate(atom_map)
        local_index[iat] = il
    end
    for l in links
        ia = local_index[l.ia]
        ib = local_index[l.ib]
        (ia > 0 || ib > 0) || continue
        push!(active_links, ActiveLink{D,T}(l, ia, ib, cl_positions[l.ia], cl_positions[l.ib]))
    end
    return active_links
end

#
# The links of the active subsystem, with the positions of the active atoms (`positions`,
# local indices) and the gradients added to the active atoms
#
function _active_link_fg!(fg::InteratomicDistanceFG{D,T}, links::Vector{ActiveLink{D,T}}, positions, sys::PackmolSystem{D,T}) where {D,T}
    fg.max_link_deviation = zero(T)
    isempty(links) && return fg
    has_pbc, unitcell, _ = _typed_unitcell(sys)
    maxdev = zero(T)
    for al in links
        l = al.link
        xa = al.ia > 0 ? positions[al.ia] : al.xa
        xb = al.ib > 0 ? positions[al.ib] : al.xb
        v = _link_vector(xa, xb, has_pbc, unitcell)
        p, g, dev = _link_penalty(v, l)
        fg.f += p
        fg.fmol[l.ma] += p
        fg.fmol[l.mb] += p
        al.ia > 0 && (fg.gxcar[al.ia] -= g)
        al.ib > 0 && (fg.gxcar[al.ib] += g)
        maxdev = max(maxdev, dev)
    end
    fg.max_link_deviation = maxdev
    return fg
end
