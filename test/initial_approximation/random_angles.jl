@testitem "_random_angles: uniform rotations" begin
    using Packmol
    using StaticArrays
    using Statistics: mean
    using Random: Xoshiro
    # The components of a rotated unit vector are uniform in [-1, 1] for uniform
    # rotations (mean absolute value 1/2). Uniform Euler angles, which this replaced,
    # gave 0.405 for the z component of the rotated z axis
    rng = Xoshiro(1)
    R = [Packmol.eulermat(Packmol._random_angles(rng, SVector{3,Float64})) for _ in 1:100_000]
    for axis in 1:3, component in 1:3
        @test abs(mean(abs(r[component, axis]) for r in R) - 0.5) < 0.01
    end
    # In two dimensions, the angles are uniform
    a = [Packmol._random_angles(rng, SVector{2,Float64}) for _ in 1:1000]
    @test all(all(0 .<= x .< 2π) for x in a)
end

@testitem "_extended_cm_bounds" setup=[TestFiles] begin
    using Packmol
    using StaticArrays
    using LinearAlgebra: norm, Diagonal
    water = TestFiles.test_file("structure_files/water.pdb")
    st = structure_type(water; number = 2, constraints = [InsideSphere([0.0, 0.0, 0.0], 20.0)])
    system = PackmolSystem([st]; output = "")
    Packmol._center_reference_coordinates!(system)
    st = system.structure_types[1]
    r = maximum(norm, st.reference_coordinates)
    lo, hi = Packmol._extended_cm_bounds(system, st, SVector(-5.0, -5.0, -5.0), SVector(5.0, 5.0, 5.0))
    @test lo ≈ SVector(-5.0, -5.0, -5.0) .- r
    @test hi ≈ SVector(5.0, 5.0, 5.0) .+ r
    # Within the periodic cell
    center = SVector(5.0, 5.0, 5.0)
    system = PackmolSystem([structure_type(water; number = 2)]; output = "",
        unitcell = Matrix{Float64}(Diagonal([10.0, 10.0, 10.0])), unitcell_center = center)
    Packmol._center_reference_coordinates!(system)
    st = system.structure_types[1]
    lo, hi = Packmol._extended_cm_bounds(system, st, SVector(0.5, 0.5, 0.5), SVector(9.5, 9.5, 9.5))
    h = Packmol._cm_half_extent(system)
    @test lo ≈ center .- h
    @test hi ≈ center .+ h
end

@testitem "initial approximation: rods in a thin spherical shell, uniform directions" begin
    # Regression test: the constraint fit of the initial approximation moves each
    # molecule to the nearest point of its region, and the random trials were drawn in
    # a box whose faces cut thin regions. Rods constrained to be radial in a thin
    # spherical shell (as lipids in a nanoliposome) ended up concentrated in the
    # directions of the corners of the box. For uniform directions, the mean of the
    # largest absolute component of a unit vector is 0.831 (0.793 for the directions
    # of uniform points of a cube)
    using Packmol
    using LinearAlgebra: norm
    using Statistics: mean
    rod = tempname() * ".pdb"
    open(rod, "w") do io
        for i in 1:6
            println(io, "ATOM  ", lpad(i, 5), "  C   ROD A   1    ", lpad("0.000", 8), lpad("0.000", 8),
                lpad(string(round(1.5 * (i - 1); digits = 3)), 8), "  1.00  0.00           C")
        end
        println(io, "END")
    end
    n = 400
    st = structure_type(rod; number = n,
        atom_constraints = [[1] => [InsideSphere([0.0, 0.0, 0.0], 20.0)], [6] => [OutsideSphere([0.0, 0.0, 0.0], 26.5)]],
    )
    stat = Float64[]
    for seed in 1:6
        system = PackmolSystem([st]; output = "", tolerance = 2.0, check = true, seed)
        redirect_stdout(devnull) do
            packmol(system)
        end
        atoms = Packmol.get_atoms(system)
        cms = [sum([a.x, a.y, a.z] for a in atoms[6(k-1)+1:6k]) / 6 for k in 1:n]
        push!(stat, mean(maximum(abs.(cm ./ norm(cm))) for cm in cms))
    end
    # 0.831 with the current code, 0.821 before (the standard error of the mean of 6
    # seeds is about 0.0025)
    @test mean(stat) > 0.826
    rm(rod; force = true)
end
