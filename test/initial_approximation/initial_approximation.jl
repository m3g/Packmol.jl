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
