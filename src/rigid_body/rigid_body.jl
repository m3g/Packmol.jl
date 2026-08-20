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
# The following function takes the vector of coordinates of a molecule, and rotates
# and translates it to put the center of mass in the origin and the principal moment
# of inertia (assuming all equal masses) along the z-axis. That will define the orientation
# of the reference coordinates.
#
function set_reference_coordinates!(x::AbstractVector{<:SVector{3,T}}; fixed::Tuple{Bool, String, Vector{T}}) where {T<:Real}
    # Move center of mass to the origin
    cm = mean(x)
    x .= x .- Ref(cm)
    # If the molecule is fixed, apply the transformation that sets its position, and return
    if first(fixed)
        # Rotate molecule according to parameters of fixed position
        rotate!(x, fixed[3][4], fixed[3][5], fixed[3][6])
        if fixed[2] == "absolute"
            x .= x + Ref(SVector(fixed[3][1:3]...))
        end
        return x
    end
    # If the molecule is not fixed, rotate it to align the principal moments of inertia
    # Inertia tensor
    I = zeros(MMatrix{3,3,T,9})
    for xi in x
        I[1,1] += xi[2]^2 + xi[3]^2
        I[2,2] += xi[1]^2 + xi[3]^2
        I[3,3] += xi[1]^2 + xi[2]^2
        I[1,2] -= xi[1] * xi[2]
        I[1,3] -= xi[1] * xi[3]
        I[2,3] -= xi[2] * xi[3]
    end
    I[2,1] = I[1,2]
    I[3,1] = I[1,3]
    I[3,2] = I[2,3]
    I = SMatrix(I)
    # Diagonalize
    evals, evecs = eigen(I)
    # Sort eigenvectors by eigenvalues
    idx = sortperm(evals, rev=true)
    evecs = evecs[:, idx]
    # Rotate
    for i in eachindex(x)
        x[i] = evecs' * x[i]
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
