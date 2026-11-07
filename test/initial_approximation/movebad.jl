@testitem "movebad!: one quota for the whole system" setup=[TestFiles] begin
    # Ten structure types of 6 molecules must not move more molecules than a
    # single type of 60 would (movefrac of the free molecules), even when every
    # molecule is bad.
    using Packmol
    using Random

    water = TestFiles.test_file("structure_files/water.pdb")
    sts = [structure_type(water; number=6, constraints=[InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])]) for _ in 1:10]
    sys = PackmolSystem(sts; output="", check_initial_placement=true, seed=1)
    Packmol.packmol(sys)
    fmol = ones(sys.nmols)
    free_mol_indices = collect(1:sys.nmols)
    mol_structure_type = Packmol._build_mol_structure_type(sys)
    rng = Random.Xoshiro(1)
    for _ in 1:20
        moved = Packmol.movebad!(sys, fmol, free_mol_indices, mol_structure_type, rng; movefrac=0.05)
        @test 1 <= length(moved) <= 3
        @test allunique(moved)
    end
end
