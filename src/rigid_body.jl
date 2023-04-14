"""
    eulermat(beta, gamma, theta, deg::String)

$(INTERNAL)

## Extended help

This routine was added because it defines the rotation in the "human" way, an is thus used
to set the position of the fixed molecules. `deg` can only be `"degree"`, in which
case the angles with be considered in degrees. If no `deg` argument
is provided, radians are used.
That means: `beta` is a counterclockwise rotation around `x` axis.
            `gamma` is a counterclockwise rotation around `y` axis.
            `theta` is a counterclockwise rotation around `z` axis.

"""
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

@testitem "eulermat" begin
    @test Packmol.eulermat(0.0, 0.0, 0.0) ≈ [1 0 0; 0 1 0; 0 0 1]
    @test Packmol.eulermat(π, 0.0, 0.0) ≈ [1 0 0; 0 -1 0; 0 0 -1]
    @test Packmol.eulermat(0.0, π, 0.0) ≈ [-1 0 0; 0 1 0; 0 0 -1]
    @test Packmol.eulermat(0.0, 0.0, π) ≈ [-1 0 0; 0 -1 0; 0 0 1]
end

"""
    move!(x::AbstractVector, newcm::AbstractVector, beta, gamma, theta)

$(INTERNAL)

## Extended help

Translates and rotates a molecule according to the desired input center of coordinates and Euler rotations modifyies the vector x.

"""
function move!(x::AbstractVector{T}, newcm::T, beta, gamma, theta) where {T<:SVector}
    cm = mean(x)
    A = eulermat(beta, gamma, theta)
    for i in eachindex(x)
        x[i] = A * (x[i] - cm) + newcm
    end
    return x
end

@testitem "move!" begin
    using StaticArrays
    x = [SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)]
    @test Packmol.move!(copy(x), SVector(0.0, 0.0, 0.0), 0.0, 0.0, 0.0) ≈
          SVector{3,Float64}[[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]
    @test Packmol.move!(copy(x), SVector(1.0, 1.0, 1.0), 0.0, 0.0, 0.0) ≈
          SVector{3,Float64}[[1.5, 1.0, 1.0], [0.5, 1.0, 1.0]]
    @test Packmol.move!(copy(x), SVector(0.0, 0.0, 0.0), π, 0.0, 0.0) ≈
          SVector{3,Float64}[[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]
    @test Packmol.move!(copy(x), SVector(0.0, 0.0, 0.0), 0.0, π, 0.0) ≈
          SVector{3,Float64}[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
    @test Packmol.move!(copy(x), SVector(0.0, 0.0, 0.0), 0.0, 0.0, π) ≈
          SVector{3,Float64}[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
end

"""
    random_move!(x_ref::AbstractVector{T}, 
                 irefatom::Int,
                 unitcell,
                 x_new::AbstractVector{T}, RNG) where {T<:SVector}

$(INTERNAL)

## Extended help

Function that generates a new random position for a molecule.
The new position is returned in `x_new`, a previously allocated array.

"""
function random_move!(
    x::AbstractVector{<:SVector{3}},
    irefatom::Int,
    unitcell,
    RNG,
)
    # To avoid boundary problems, the center of coordinates are generated in a 
    # much larger region, and wrapped aftwerwards
    scale = 100.0

    # Generate random coordinates for the center of mass
    cmin, cmax = PeriodicSystems.get_computing_box(system)
    newcm = SVector{3}(scale * (cmin[i] + random(RNG, Float64) * (cmax[i] - cmin[i])) for i in 1:3)

    # Generate random rotation angles 
    beta = 2π * random(RNG, Float64)
    gamma = 2π * random(RNG, Float64)
    theta = 2π * random(RNG, Float64)

    # Take care that this molecule is not split by periodic boundary conditions, by
    # wrapping its coordinates around its reference atom
    for iat in eachindex(x)
        x[iat] = CellListMap.wrap_relative_to(x[iat], x[irefatom], unitcell)
    end

    # Move molecule to new position
    move!(x, newcm, beta, gamma, theta)

    return x
end

@testitem "random_move!" begin
    using Packmol
    using StaticArrays
    using LinearAlgebra: norm
    using CellListMap.PeriodicSystems
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

    RNG = Packmol.init_random(Options())
    # Orthorhombic cell
    x = [-1.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    system = PeriodicSystem(positions=x, cutoff=0.1, unitcell=SVector(10.0, 10.0, 10.0), output=0.0)
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [-9.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [4.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))

    # Triclinic cell
    x = [-1.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    system = PeriodicSystem(positions=x, cutoff=0.1, unitcell=@SMatrix[10.0 5.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0], output=0.0)
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [-9.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [4.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))

end