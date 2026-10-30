#
# Rotation constraints (`constrain_rotation`): the angle between an axis (x, y or z)
# and a vector of the molecule must be within `lower` and `upper`. The vector joins the
# centers of two sets of atoms of the molecule (`from` and `to`, indices of the atoms of
# the structure), or, when no sets are given, it is the longest principal axis of the
# molecule (z of its reference orientation, see `set_reference_coordinates!`).
#
# The vector is a fixed linear combination of the atom positions of the molecule,
# v = Σ wᵢ xᵢ with Σ wᵢ = 0: wᵢ = 1/n_to on the `to` atoms and -1/n_from on the `from`
# atoms; for the principal axis, wᵢ = L zᵢ / Σ zⱼ² (z of the reference coordinates, L
# the length of the molecule along z), for which v = L R ẑ exactly for a rigid molecule. The penalty, `weight * (L * Δθ)^2`, with
# Δθ the violation of the angle and L = |v|, is the squared displacement of the tip of
# the vector, comparable to the penalties of the other constraints. Its gradient is
# applied to the atoms through the weights, so it goes through the same chain rule.
# L is constant for a rigid molecule, so the gradient with respect to v is the one at
# constant |v| (perpendicular to v).
#
weight_default[:rotation] = 5.0

mutable struct RotationConstraint{D,T}
    axis::SVector{D,T}
    lower::T # radians, in [0, π]
    upper::T
    from::Vector{Int}
    to::Vector{Int}
    weight::T
    # Set from the reference coordinates by `_set_rotation_constraint!`
    indices::Vector{Int}
    weights::Vector{T}
    reference_vector::SVector{D,T}
end

function RotationConstraint{D,T}(
    axis::Symbol, center_deg, halfwidth_deg;
    from::AbstractVector{<:Integer}=Int[], to::AbstractVector{<:Integer}=Int[],
    weight=weight_default[:rotation],
) where {D,T}
    iaxis = axis == :x ? 1 : axis == :y ? 2 : axis == :z ? 3 : 0
    (1 <= iaxis <= D) || throw(ArgumentError("""\n
        constrain_rotation axis must be $(D == 3 ? "x, y or z" : "x or y"), got $axis.

    """))
    isempty(from) == isempty(to) || throw(ArgumentError("""\n
        constrain_rotation: both sets of atoms (`atoms` and `to`) must be given, or none.

    """))
    c = T(center_deg) * T(π) / 180
    hw = abs(T(halfwidth_deg)) * T(π) / 180
    return RotationConstraint{D,T}(
        SVector{D,T}(ntuple(i -> i == iaxis ? one(T) : zero(T), D)),
        clamp(c - hw, zero(T), T(π)), clamp(c + hw, zero(T), T(π)),
        collect(Int, from), collect(Int, to), T(weight),
        Int[], T[], zero(SVector{D,T}),
    )
end

#
# Weights of the atoms in the vector of the constraint, from the reference coordinates
# of the molecule (centered, and in their principal axes, for free molecules).
#
function _set_rotation_constraint!(rc::RotationConstraint{D,T}, ref::AbstractVector{SVector{D,T}}) where {D,T}
    natoms = length(ref)
    if isempty(rc.from)
        s = sum(r -> r[D]^2, ref)
        s > sqrt(eps(T)) || throw(ArgumentError("""\n
            constrain_rotation without sets of atoms needs a molecule with a longest axis:
            give the sets of atoms that define the vector.

        """))
        L = maximum(r -> r[D], ref) - minimum(r -> r[D], ref)
        rc.indices = collect(1:natoms)
        rc.weights = [L * r[D] / s for r in ref]
    else
        for i in Iterators.flatten((rc.from, rc.to))
            (1 <= i <= natoms) || throw(ArgumentError("""\n
                constrain_rotation: atom index $i out of range (the molecule has $natoms atoms).

            """))
        end
        rc.indices = vcat(rc.from, rc.to)
        rc.weights = vcat(fill(-one(T) / length(rc.from), length(rc.from)), fill(one(T) / length(rc.to), length(rc.to)))
    end
    v = zero(SVector{D,T})
    for (i, w) in zip(rc.indices, rc.weights)
        v += w * ref[i]
    end
    norm(v) > sqrt(eps(T)) || throw(ArgumentError("""\n
        constrain_rotation: the two sets of atoms have the same center.

    """))
    rc.reference_vector = v
    return rc
end

#
# Penalty of the rotation constraint for the vector `v` (`R * reference_vector`), and its
# gradient with respect to `v`.
#
function rotation_penalty(rc::RotationConstraint{D,T}, v::SVector{D,T}) where {D,T}
    L = norm(v)
    θ = acos(clamp(dot(v, rc.axis) / L, -one(T), one(T)))
    Δ = θ < rc.lower ? θ - rc.lower : θ > rc.upper ? θ - rc.upper : zero(T)
    return rc.weight * (L * Δ)^2
end

function rotation_gradient(rc::RotationConstraint{D,T}, v::SVector{D,T}) where {D,T}
    L = norm(v)
    cosθ = clamp(dot(v, rc.axis) / L, -one(T), one(T))
    θ = acos(cosθ)
    Δ = θ < rc.lower ? θ - rc.lower : θ > rc.upper ? θ - rc.upper : zero(T)
    Δ == zero(T) && return zero(SVector{D,T})
    sinθ = sqrt(max(zero(T), one(T) - cosθ^2))
    # dθ/dv = -(axis - cos θ v̂) / (L sin θ); undefined (zero) at θ = 0 or π
    sinθ < sqrt(eps(T)) && return zero(SVector{D,T})
    dθdv = -(rc.axis - cosθ * v / L) / (L * sinθ)
    return 2 * rc.weight * L^2 * Δ * dθdv
end

# Total rotation penalty of a molecule of the structure type `st` with rotation `R`.
function _rotation_penalty(st, R::SMatrix{D,D,T}) where {D,T}
    p = zero(T)
    for rc in st.rotation_constraints
        p += rotation_penalty(rc, R * rc.reference_vector)
    end
    return p
end

#
# Adds the rotation penalties of the molecule whose first atom is `iat_first`, with
# rotation `R`, to the Cartesian gradient `gxcar`, and returns the penalty.
#
function _rotation_fg!(gxcar::AbstractVector{SVector{D,T}}, st, R::SMatrix{D,D,T}, iat_first::Int) where {D,T}
    p = zero(T)
    for rc in st.rotation_constraints
        v = R * rc.reference_vector
        p += rotation_penalty(rc, v)
        g = rotation_gradient(rc, v)
        g == zero(SVector{D,T}) && continue
        for (i, w) in zip(rc.indices, rc.weights)
            gxcar[iat_first+i-1] += w * g
        end
    end
    return p
end
