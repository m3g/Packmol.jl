@testitem "eulermat" begin
    using StaticArrays
    @test Packmol.eulermat(0.0, 0.0, 0.0) ≈ [1 0 0; 0 1 0; 0 0 1]
    @test Packmol.eulermat(π, 0.0, 0.0) ≈ [1 0 0; 0 -1 0; 0 0 -1]
    @test Packmol.eulermat(0.0, π, 0.0) ≈ [-1 0 0; 0 1 0; 0 0 -1]
    @test Packmol.eulermat(0.0, 0.0, π) ≈ [-1 0 0; 0 -1 0; 0 0 1]
    @test Packmol.eulermat(SVector(0.0, 0.0, π)) ≈ [-1 0 0; 0 -1 0; 0 0 1]
end

#=
    move!(x::AbstractVector, newcm::AbstractVector, beta, gamma, theta)

$(INTERNAL)

## Extended help

Translates and rotates a molecule according to the desired input center of coordinates and Euler rotations modifyies the vector x.

=#

@testitem "move!" setup=[RigidBody] begin
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
    @test Packmol.move!(copy(x), Packmol.MoleculePosition(SVector(0.0, 0.0, 0.0), SVector(0.0, π, 0.0))) ≈
          SVector{3,Float64}[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
    @test Packmol.move!(copy(x), Packmol.MoleculePosition(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, π))) ≈
          SVector{3,Float64}[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
end

#=
    random_move!(x_ref::AbstractVector{T}, 
                 irefatom::Int,
                 unitcell,
                 x_new::AbstractVector{T}, RNG) where {T<:SVector}

$(INTERNAL)

## Extended help

Function that generates a new random position for a molecule.
The new position is returned in `x_new`, a previously allocated array.

=#

@testitem "random_move!" setup=[RigidBody] begin
    # Orthorhombic cell
    x = [-1.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    system = ParticleSystem(positions=x, cutoff=0.1, unitcell=SVector(10.0, 10.0, 10.0), output=0.0)
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [-9.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [4.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    # Triclinic cell
    x = [-1.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    system = ParticleSystem(positions=x, cutoff=0.1, unitcell=@SMatrix[10.0 5.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0], output=0.0)
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [-9.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
    system.xpositions .= [4.0 .+ 2 * rand(SVector{3,Float64}) for _ = 1:5]
    @test check_internal_distances(x, Packmol.random_move!(copy(x), 1, system, RNG))
end

@testitem "set_reference_coordinates!" setup=[RigidBody] begin
    using PDBTools
    using StaticArrays
    using BenchmarkTools
    using Statistics: mean
    # A molecule not `fixed` in the input file: `set_reference_coordinates!`
    # still requires the keyword, but never reads past `first(fixed)` in that case.
    not_fixed = (false, "", Float64[])
    x = [ SVector(0.0, 0.0, 0.0), SVector(√3/3, √3/3, √3/3) ]
    Packmol.set_reference_coordinates!(x; fixed=not_fixed)
    @test x[1] ≈ SVector(0., 0., -0.5)
    @test x[2] ≈ SVector(0., 0., 0.5)
    x = [ SVector(√3/3, √3/3, √3/3), SVector(0.0, 0.0, 0.0) ]
    Packmol.set_reference_coordinates!(x; fixed=not_fixed)
    @test x[1] ≈ SVector(0.0, 0.0, 0.5)
    @test x[2] ≈ SVector(0.0, 0.0, -0.5)
    water = coor(read_pdb(joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")))
    water_save = copy(water)
    not_fixed_water = (false, "", zeros(eltype(eltype(water)), 0))
    a = @ballocated Packmol.set_reference_coordinates!($water; fixed=$not_fixed_water) samples=1 evals=1
    @test a == 0
    # water coordinates are Float32 (as read from the PDB file), so the
    # residual from centering is limited by Float32 precision, not 1e-10.
    @test mean(water) ≈ SVector(0.0, 0.0, 0.0) atol = 1e-6
    @test all(isapprox(v[1],0.0,atol=1e-10) for v in water)
    @test check_internal_distances(water, water_save)
end
