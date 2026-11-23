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

@testitem "movebad!: overlap of the trial positions with the current atoms" setup=[TestFiles] begin
    # The score of the trial positions of movebad: the overlap of a molecule, at a given
    # position, with the atoms of the others at their current positions
    using Packmol
    using Random
    using StaticArrays

    water = TestFiles.test_file("structure_files/water.pdb")
    st = structure_type(water; number=20, constraints=[InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])])
    sys = PackmolSystem([st]; output="", check_initial_placement=true, seed=1)
    Packmol.packmol(sys)
    atom_positions = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    current = Packmol._current_atoms(sys, atom_positions, sys.structure_types)
    ref = sys.structure_types[1].reference_coordinates
    mp1 = sys.molecule_positions[1]
    mp2 = sys.molecule_positions[2]
    # molecule 1 on top of molecule 2: overlaps it
    @test Packmol._current_overlap(current, mp2, ref, 1, 0) > 0
    # skipping molecule 2 itself: only the overlaps with the others
    @test Packmol._current_overlap(current, mp2, ref, 2, 0) < Packmol._current_overlap(current, mp2, ref, 1, 0)
    # far from all the molecules: no overlap
    far = Packmol.MoleculePosition(mp1.cm .+ 1000.0, mp1.angles)
    @test Packmol._current_overlap(current, far, ref, 1, 0) == 0
    # movebad with the score: the moved molecules satisfy their constraints
    fmol = ones(sys.nmols)
    mol_structure_type = Packmol._build_mol_structure_type(sys)
    cm_min, cm_max = Packmol.compute_cm_bounds(sys)
    current = Packmol._current_atoms(sys, atom_positions, sys.structure_types)
    moved = Packmol.movebad!(sys, fmol, collect(1:sys.nmols), mol_structure_type, Random.Xoshiro(1);
        movefrac=0.1, cm_lo_type=cm_min, cm_hi_type=cm_max, current, ntrials=Packmol._movebad_trials)
    @test !isempty(moved)
    for imol in moved
        @test all(0.0 .<= sys.molecule_positions[imol].cm .<= 20.0)
    end
end
