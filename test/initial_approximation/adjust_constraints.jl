@testitem "adjust_constraints: state after a failed adjustment" setup = [TestFiles] begin
    # Regression test: when the constraint adjustment of the initial approximation does
    # not converge, its last loop re-placed the bad molecules at random positions that
    # were not fitted again; and a structure type of a single molecule (degenerate bounds
    # of its center of mass) was re-placed in the box of ±sidemax around the origin, far
    # from a system that is not centered there. The packing then started from atoms
    # 2000 Å out of their sphere (the long RNAs of a cell).
    using Packmol
    using LinearAlgebra: norm
    water = TestFiles.test_file("structure_files/water.pdb")
    c = [3000.0, 3000.0, 3000.0]
    types = [
        structure_type(water; number = 50, constraints = [InsideSphere(c, 15.0)]),
        # impossible: inside a sphere of 2 Å and outside one of 4 Å, with the same center
        structure_type(water; number = 1, constraints = [InsideSphere(c, 2.0), OutsideSphere(c, 4.0)]),
    ]
    system = PackmolSystem(types; output = "", tolerance = 2.0, seed = 1, nloop = 5,
        check_initial_placement = true, continue_failed_initial_packing = true)
    redirect_stdout(devnull) do
        packmol(system)
    end
    @test norm(system.molecule_positions[end].cm - c) < 10.0
    @test count(mp -> norm(mp.cm - c) > 20.0, system.molecule_positions) == 0
    # The union of the valid bounds
    lo, hi = Packmol._union_of_valid_bounds(
        [Packmol.SVector(0.0, 0.0, 0.0), Packmol.SVector(1.0, 1.0, 1.0), Packmol.SVector(-5.0, 2.0, 2.0)],
        [Packmol.SVector(1.0, 2.0, 3.0), Packmol.SVector(1.0, 1.0, 1.0), Packmol.SVector(0.0, 4.0, 2.5)],
    )
    @test lo == Packmol.SVector(-5.0, 0.0, 0.0) && hi == Packmol.SVector(1.0, 4.0, 3.0)
    @test Packmol._union_of_valid_bounds(nothing, nothing) == (nothing, nothing)
end
