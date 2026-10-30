#=
    eulermat(beta, gamma, theta, deg::String)

$(INTERNAL)

## Extended help

This routine returns a rotation matrix that rotates a vector by `beta`, `gamma`, and `theta` angles.

It defines the rotation in the "human" way, an is thus used to set the position of the fixed molecules. `deg` can only be `"degree"`, in which
case the angles with be considered in degrees. If no `deg` argument
is provided, radians are used.

That means: `beta` is a counterclockwise rotation around `x` axis.
            `gamma` is a counterclockwise rotation around `y` axis.
            `theta` is a counterclockwise rotation around `z` axis.

=#
function eulermat(beta, gamma, theta, deg::String)
    if deg != "degree"
        error("ERROR: to use radians just omit the last parameter")
    end
    beta = beta * π / 180
    gamma = gamma * π / 180
    theta = theta * π / 180
    return eulermat(beta, gamma, theta)
end

function eulermat(beta, gamma, theta)
    c1 = cos(beta)
    s1 = sin(beta)
    c2 = cos(gamma)
    s2 = sin(gamma)
    c3 = cos(theta)
    s3 = sin(theta)
    #! format: off
    @SMatrix [
               c2*c3          -c2*s3         s2
        (c1*s3+c3*s1*s2) (c1*c3-s1*s2*s3) -c2*s1
        (s1*s3-c1*c3*s2) (c1*s2*s3+c3*s1)  c1*c2
    ]
    #! format: on
end

eulermat(angs::SVector{3}) = eulermat(angs[1], angs[2], angs[3])

function move!(x::AbstractVector{T}, newcm::T, beta, gamma, theta) where {T<:SVector}
    cm = mean(x)
    x .= x .- Ref(cm)
    rotate!(x, beta, gamma, theta)
    x .= x .+ Ref(newcm)
    return x
end
move!(x::AbstractVector{<:SVector{3,T}}, newcm::SVector{3,T}, angles::SVector{3,T}) where {T} =
    move!(x, newcm, angles[1], angles[2], angles[3])
move!(x::AbstractVector{<:SVector{D,T}}, pos::MoleculePosition{D,T}) where {D,T} =
    move!(x, pos.cm, pos.angles)

function rotate!(x::AbstractVector{T}, beta, gamma, theta) where {T<:SVector}
    A = eulermat(beta, gamma, theta)
    for i in eachindex(x)
        x[i] = A * x[i]
    end
    return x
end
rotate!(x::AbstractVector{<:SVector{D,T}}, pos::MoleculePosition{D,T}) where {D,T} =
    rotate!(x, pos.angles)
rotate!(x::AbstractVector{<:SVector{3,T}}, angles::SVector{3,T}) where {T} =
    rotate!(x, angles[1], angles[2], angles[3])

#
# Optimal rigid-body rotation (Kabsch algorithm) aligning `reference` onto
# `observed`, both already centered at their own centroid: the rotation
# matrix R minimizing `sum(norm(R*reference[i] - observed[i])^2 for i)`.
# `reference`/`observed` must correspond atom-for-atom (same length, same
# order). Used to recover a molecule's (cm, angles) from restart coordinates
# — see `_align_molecule` in write_output.jl.
#
function kabsch_rotation(reference::AbstractVector{<:SVector{3,T}}, observed::AbstractVector{<:SVector{3,T}}) where {T}
    length(reference) == length(observed) ||
        throw(ArgumentError("reference and observed must have the same number of points, got $(length(reference)) and $(length(observed))"))
    H = sum(r * o' for (r, o) in zip(reference, observed))
    F = svd(Matrix(H))
    d = sign(det(F.V * F.U'))
    D = Diagonal(SVector(one(T), one(T), T(d)))
    return SMatrix{3,3,T}(F.V * D * F.U')
end

#
# Inverse of `eulermat`: recovers (beta, gamma, theta) — in eulermat's own
# convention (x, y, z rotation, respectively) — from a rotation matrix, up
# to the usual gimbal-lock ambiguity at gamma = ±π/2 (there beta and theta
# are not individually identifiable from R alone; beta is arbitrarily fixed
# to 0 and theta absorbs the whole remaining rotation).
#
function euler_angles(R::AbstractMatrix{T}) where {T}
    s2 = clamp(R[1, 3], -one(T), one(T))
    gamma = asin(s2)
    c2 = sqrt(one(T) - s2^2)
    # A generous threshold, not just "not exactly zero": `R[2,3]`/`R[3,3]`
    # (used in the well-conditioned branch below) both scale with c2, so
    # their ratio becomes numerically unstable well before c2 actually
    # reaches zero — particularly for an `R` that, unlike a literal
    # `eulermat(...)` call, is itself already the result of a noisy
    # computation (e.g. `kabsch_rotation`'s SVD).
    if c2 > sqrt(sqrt(eps(T)))
        beta = atan(-R[2, 3], R[3, 3])
        theta = atan(-R[1, 2], R[1, 1])
    else
        beta = zero(T)
        theta = atan(R[2, 1], R[2, 2])
    end
    return SVector{3,T}(beta, gamma, theta)
end

function random_move!(
    x::AbstractVector{<:SVector{3}},
    irefatom::Int,
    system,
    RNG,
)
    # To avoid boundary problems, the center of coordinates are generated in a 
    # much larger region, and wrapped aftwerwards
    scale = 100.0
    # Generate random coordinates for the center of mass
    cmin, cmax = CellListMap.get_computing_box(system)
    newcm = SVector{3}(scale * (cmin[i] + rand(RNG, Float64) * (cmax[i] - cmin[i])) for i in 1:3)
    # Generate random rotation angles 
    beta = 2π * rand(RNG, Float64)
    gamma = 2π * rand(RNG, Float64)
    theta = 2π * rand(RNG, Float64)
    # Take care that this molecule is not split by periodic boundary conditions, by
    # wrapping its coordinates around its reference atom
    for iat in eachindex(x)
        x[iat] = CellListMap.wrap_relative_to(x[iat], x[irefatom], system.unitcell)
    end
    # Move molecule to new position
    move!(x, newcm, beta, gamma, theta)
    return x
end

#
# Center the coordinates of a molecule at the origin and rotate them to the principal
# axes of inertia of the molecule (with equal masses): the axis of the smallest moment
# of inertia (along which the molecule is longest) along the last axis (z in 3D), and
# that of the largest moment along x. This is the reference orientation of the free
# molecules, the orientation of zero Euler angles. Each axis is chosen with its largest
# component positive, and the rotation is proper, so that aligning coordinates that are
# already aligned does not change them.
#
_largest_component_sign(v::SVector{D,T}) where {D,T} = v[argmax(abs.(v))] < 0 ? -one(T) : one(T)

function set_reference_coordinates!(x::AbstractVector{<:SVector{D,T}}) where {D,T<:Real}
    x .= x .- Ref(mean(x))
    C = zero(SMatrix{D,D,T})
    for xi in x
        C += xi * xi'
    end
    # Eigenvalues in increasing order (the spread along each axis): the longest axis last
    evecs = eigen(Symmetric(C)).vectors
    signs = SVector{D,T}(ntuple(i -> _largest_component_sign(evecs[:, i]), Val(D)))
    R = evecs .* signs'
    if det(R) < 0
        R = R .* SVector{D,T}(ntuple(i -> i == 1 ? -one(T) : one(T), Val(D)))'
    end
    for i in eachindex(x)
        x[i] = R' * x[i]
    end
    return x
end

@testsnippet RigidBody begin
    using Packmol
    using CellListMap
    using StaticArrays
    using LinearAlgebra: norm
    import Random
    RNG = Random.Xoshiro()
    # Function that checks that the internal distances of a molecule are preserved
    function check_internal_distances(x, y)
        for i = firstindex(x):lastindex(x)-1
            for j = i+1:lastindex(x)
                x_ij = norm(x[i] - x[j])
                y_ij = norm(y[i] - y[j])
                if !isapprox(x_ij, y_ij)
                    return false, x_ij, y_ij
                end
            end
        end
        return true
    end
end
