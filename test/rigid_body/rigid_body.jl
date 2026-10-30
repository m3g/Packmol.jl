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

@testitem "kabsch_rotation / euler_angles" setup=[RigidBody] begin
    using StaticArrays
    using LinearAlgebra: norm
    import Random

    # A generic (5-point, non-degenerate) centered point cloud, rotated by a
    # known R, must have that R (and its Euler-angle decomposition) exactly
    # recoverable — this is what `_align_molecule` relies on for a PDB-based
    # restart_from.
    rng = Random.Xoshiro(1)
    n = 5
    pts = [SVector{3,Float64}(randn(rng, 3)) for _ in 1:n]
    pts .-= Ref(sum(pts) / n)

    for _ in 1:200
        beta, gamma, theta = 2π .* rand(rng, 3) .- π
        R_true = Packmol.eulermat(beta, gamma, theta)
        observed = [R_true * p for p in pts]

        R_rec = Packmol.kabsch_rotation(pts, observed)
        @test maximum(norm(R_rec * pts[i] - observed[i]) for i in 1:n) < 1e-8

        angles = Packmol.euler_angles(R_rec)
        R_from_angles = Packmol.eulermat(angles...)
        @test maximum(norm(R_from_angles * pts[i] - observed[i]) for i in 1:n) < 1e-3
    end

    # Gimbal lock (gamma = ±π/2): beta and theta aren't individually
    # identifiable from R alone, but the reconstructed rotation must still
    # match — see `euler_angles`'s docstring.
    for (beta, gamma, theta) in ((0.3, π / 2, 0.7), (1.1, -π / 2, -0.4), (2.0, -π / 2, -1.5))
        R_true = Packmol.eulermat(beta, gamma, theta)
        observed = [R_true * p for p in pts]
        R_rec = Packmol.kabsch_rotation(pts, observed)
        angles = Packmol.euler_angles(R_rec)
        R_from_angles = Packmol.eulermat(angles...)
        @test maximum(norm(R_from_angles * pts[i] - observed[i]) for i in 1:n) < 1e-6
    end

    @test_throws ArgumentError Packmol.kabsch_rotation(pts, pts[1:end-1])
end

@testitem "set_reference_coordinates!" setup=[RigidBody, TestFiles] begin
    using PDBToolsBase: read_pdb, positions
    using StaticArrays
    using BenchmarkTools
    using Statistics: mean
    using LinearAlgebra: det
    # Centered, with the longest axis along z
    x = [ SVector(0.0, 0.0, 0.0), SVector(√3/3, √3/3, √3/3) ]
    Packmol.set_reference_coordinates!(x)
    @test x[1] ≈ SVector(0., 0., -0.5) atol = 1e-12
    @test x[2] ≈ SVector(0., 0., 0.5) atol = 1e-12
    # Principal axes: a rotated and translated ellipsoidal cloud of points ends up with a
    # diagonal, increasing covariance; with a proper rotation (internal distances and
    # handedness kept); aligning again does not change the coordinates
    y = [SVector(cos(i), 3 * sin(2i), 0.3 * cos(5i)) for i in 1:20]
    y = [Packmol.eulermat(SVector(0.3, 1.2, -0.7)) * v .+ 4 for v in y]
    y0 = copy(y)
    Packmol.set_reference_coordinates!(y)
    C = sum(v * v' for v in y)
    @test all(abs(C[i, j]) < 1e-10 for i in 1:3, j in 1:3 if i != j)
    @test C[1, 1] < C[2, 2] < C[3, 3]
    @test check_internal_distances(y, y0)
    triple(v) = det(hcat(v[2] - v[1], v[3] - v[1], v[4] - v[1]))
    @test triple(y) ≈ triple(y0)
    y1 = copy(y)
    Packmol.set_reference_coordinates!(y)
    @test maximum(norm.(y .- y1)) < 1e-12
    # 2D
    z = [SVector(3 * cos(i), sin(3i)) for i in 1:10]
    Packmol.set_reference_coordinates!(z)
    C = sum(v * v' for v in z)
    @test abs(C[1, 2]) < 1e-10 && C[1, 1] < C[2, 2]
    water = positions(read_pdb(TestFiles.test_file("structure_files/water.pdb")))
    water_save = copy(water)
    a = @ballocated Packmol.set_reference_coordinates!($water) samples=1 evals=1
    @test a == 0
    # water coordinates are Float32 (as read from the PDB file), so the
    # residual from centering is limited by Float32 precision, not 1e-10.
    @test mean(water) ≈ SVector(0.0, 0.0, 0.0) atol = 1e-6
    # water is planar: the axis normal to its plane (the largest moment of inertia) is x
    @test all(isapprox(v[1], 0.0, atol=1e-5) for v in water)
    @test check_internal_distances(water, water_save)
end
