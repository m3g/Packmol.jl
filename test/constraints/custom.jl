@testitem "CustomConstraint penalty and gradient" begin
    using StaticArrays, LinearAlgebra, ForwardDiff
    # A sphere of radius 5: penalty 5 (|x| - 5)^2 outside it
    f(x) = norm(x) - 5.0
    penalty(x) = f(x) > 0 ? 5.0 * f(x)^2 : 0.0
    gradient(x) = f(x) > 0 ? 10.0 * f(x) * x / norm(x) : zero(x)
    for c in (
        CustomConstraint(f),                                            # finite differences
        CustomConstraint(f; gradient = x -> ForwardDiff.gradient(f, x)), # automatic differentiation
        CustomConstraint(f; gradient = x -> x / norm(x)),               # analytical
    )
        for x in (SVector(1.0, 2.0, 3.0), SVector(4.0, -3.0, 6.0), SVector(-10.0, 0.5, 2.0))
            @test Packmol.constraint_penalty(c, x) ≈ penalty(x)
            @test isapprox(Packmol.constraint_gradient(c, x), gradient(x); rtol = 1e-6)
        end
        @test Packmol.constraint_gradient(c, SVector(1.0, 1.0, 1.0)) == zeros(3)
    end
    # The weight, and the type of the result of the user function
    c = CustomConstraint(x -> Float32(x[3]); weight = 2.0)
    @test Packmol.constraint_penalty(c, SVector(0.0, 0.0, 3.0)) ≈ 2.0 * 9.0
    @test occursin("CustomConstraint", sprint(show, MIME"text/plain"(), c))
end

@testitem "CustomConstraint packing" setup = [TestFiles] begin
    using StaticArrays, LinearAlgebra
    water = TestFiles.test_file("structure_files/water.pdb")
    # Inside a torus of radii 12 and 5 Å around the z axis, and the oxygen atoms below z = 2
    f(x) = hypot(hypot(x[1], x[2]) - 12.0, x[3]) - 5.0
    st = structure_type(water; number = 60,
        constraints = [CustomConstraint(f)],
        atom_constraints = [[1] => CustomConstraint(x -> x[3] - 2.0)],
    )
    sys = PackmolSystem([st]; tolerance = 2.0, output = tempname() * ".pdb")
    redirect_stdout(devnull) do
        packmol(sys)
    end
    @test sys.status == :packing_ready
    atoms = get_atoms(sys)
    @test maximum(f(SVector(a.x, a.y, a.z)) for a in atoms) < 0.1
    @test maximum(a.z for a in atoms if a.name == "OH2" || startswith(a.name, "O")) < 2.1
end
