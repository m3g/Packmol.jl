@testitem "tolerance_precision: 5% of twice the smallest radius" setup=[TestFiles] begin
    file = TestFiles.test_file("structure_files/water.pdb")
    box = InsideBox([0., 0., 0.], [20., 20., 20.])
    # default radii (half the tolerance): 5% of the tolerance
    st = structure_type(file; number = 10, constraints = [box])
    sys = PackmolSystem([st]; output = tempname() * ".pdb", tolerance = 2.0)
    @test Packmol._tolerance_precision(sys) ≈ 0.1
    sys = PackmolSystem([structure_type(file; number = 10, constraints = [box], tolerance = 3.0)];
        output = tempname() * ".pdb", tolerance = 3.0)
    @test Packmol._tolerance_precision(sys) ≈ 0.15
    # the smallest radius of all atoms
    small = structure_type(file; number = 10, constraints = [box], radius = 0.6)
    large = structure_type(file; number = 1, fixed = ([0., 0., 0.], [0., 0., 0.]), radius = 1.5)
    sys = PackmolSystem([small, large]; output = tempname() * ".pdb", tolerance = 2.0)
    @test Packmol._tolerance_precision(sys) ≈ 0.06
    sys = PackmolSystem([large, structure_type(file; number = 10, constraints = [box], radius = 1.5)];
        output = tempname() * ".pdb", tolerance = 2.0)
    @test Packmol._tolerance_precision(sys) ≈ 0.15
    # set explicitly
    sys.tolerance_precision = 0.02
    @test Packmol._tolerance_precision(sys) == 0.02
    # changing the tolerance changes the default radii, and so the precision
    sys = PackmolSystem([st]; output = tempname() * ".pdb", tolerance = 2.0)
    Packmol._set_tolerance!(sys, 4.0)
    @test Packmol._tolerance_precision(sys) ≈ 0.2
end
