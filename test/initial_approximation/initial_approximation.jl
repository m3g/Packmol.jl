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
        system = PackmolSystem([st]; output = "", tolerance = 2.0, check_initial_placement = true, seed)
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

@testitem "initial approximation: a Gaussian drop inside a cylinder fills the dome" begin
    # Regression test: the region below a radial Gaussian, above a plane and inside a
    # cylinder (the water drop tutorial) was not used for the direct sampling of the
    # initial placement (cylinders were left out), so the molecules started in the large
    # box and the constraint fit pulled them to the nearest points of the region: a flat
    # layer on the plane, never higher than 9 Å, under a Gaussian 27 Å high at the axis.
    using Packmol
    point = tempname() * ".pdb"
    open(point, "w") do io
        println(io, "ATOM      1  O   WAT A   1       0.000   0.000   0.000  1.00  0.00           O")
        println(io, "END")
    end
    st = structure_type(point; number = 1000, constraints = [
        BelowRadialGaussian([0., 0., 1.], [0., 0., 0.], 2., 25., 15.),
        AbovePlane([0., 0., 1.], 2.),
        InsideCylinder([0., 0., 0.], [0., 0., 1.], 35., 40.),
    ])
    system = PackmolSystem([st]; output = "", tolerance = 2.0, check_initial_placement = true)
    redirect_stdout(devnull) do
        packmol(system)
    end
    atoms = Packmol.get_atoms(system)
    near_axis = [a.z for a in atoms if hypot(a.x, a.y) < 8]
    @test !isempty(near_axis)
    @test maximum(near_axis) > 18
    # and all within the cylinder
    @test all(hypot(a.x, a.y) < 35.1 for a in atoms)
    rm(point; force = true)
end
