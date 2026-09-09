@testitem "water box packing" begin
    using Packmol
    using PDBTools: read_pdb
    using StaticArrays
    using LinearAlgebra: norm

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 100
    @test length(sys.atoms) == 300
    @test sys.tolerance == 2.0
    @test sys.status == :not_packed

    result = Packmol.packmol(sys; nloop=200, maxit=20, iprint=10)
    @test result === sys
    @test sys.status == :packing_ready

    # Compute final atom positions
    atom_positions = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    Packmol.compute_atom_positions!(atom_positions, sys.molecule_positions, sys)

    # Check that all atoms satisfy the box constraint (penalty ≈ 0)
    st = sys.structure_types[1]
    c = st.constraints[1]  # InsideBox
    for pos in atom_positions
        @test Packmol.constraint_penalty(c, pos) ≈ 0.0 atol = 0.1
    end

    # Check minimum inter-molecular distance is close to tolerance
    natoms_per_mol = st.natoms
    dmin = let dmin = Inf
        for i in 1:length(atom_positions)
            mol_i = (i - 1) ÷ natoms_per_mol + 1
            for j in i+1:length(atom_positions)
                mol_j = (j - 1) ÷ natoms_per_mol + 1
                if mol_i != mol_j
                    d = norm(atom_positions[i] - atom_positions[j])
                    dmin = min(dmin, d)
                end
            end
        end
        dmin
    end
    @test dmin >= sys.tolerance - 0.1

    # Write output and verify
    output_file = Packmol.write_output(sys)
    @test isfile(output_file)
    output_atoms = read_pdb(output_file)
    @test length(output_atoms) == 300
    rm(output_file; force=true)
end

@testitem "constrain_rotation" begin
    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    water = structure_type(water_pdb; number=5,
        constraints=[InsideBox([0, 0, 0], [28, 28, 28])],
        constrain_rotation=Dict(:x => (0.0, 5.0), :y => (0.0, 5.0), :z => (0.0, 5.0)),
    )
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0),
        unitcell_center=zeros(3),
    )
    result = packmol(sys; nloop=20, maxit=100, iprint=1000)
    @test result === sys
    @test sys.status == :packing_ready
    for mp in sys.molecule_positions
        for a in mp.angles
            @test -5.0 * pi / 180 <= a <= 5.0 * pi / 180
        end
    end
    rm(sys.output_file; force=true)
end

@testitem "restart_from / restart_to" begin
    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    urea_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "urea.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    restart_all = tempname()
    restart_water = tempname()

    # First run: converge normally, writing whole-system and per-type restart files
    water = structure_type(water_pdb; number=4,
        constraints=[InsideBox([0, 0, 0], [28, 28, 28])], restart_to=restart_water)
    urea = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys1 = PackmolSystem([water, urea]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_to=restart_all,
    )
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=1).status == :packing_ready
    @test isfile(restart_all)
    @test isfile(restart_water)
    saved_water_cm = sys1.molecule_positions[1].cm
    rm(sys1.output_file; force=true)

    # Second run: whole-system restart_from should skip initial placement
    # entirely and reproduce the exact same (already-converged) positions —
    # converging immediately, in the very first loop.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    urea2 = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2, urea2]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=2).status == :packing_ready
    @test sys2.molecule_positions[1].cm ≈ saved_water_cm
    rm(sys2.output_file; force=true)

    # Third run: per-structure-type restart_from applies only to that type,
    # and (regression test) must NOT leak into the whole-system setting —
    # the other structure type still gets a normal (non-restarted) placement.
    water3 = structure_type(water_pdb; number=4,
        constraints=[InsideBox([0, 0, 0], [28, 28, 28])], restart_from=restart_water)
    urea3 = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([water3, urea3]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3),
    )
    @test isnothing(sys3.restart_from)
    @test sys3.structure_types[1].restart_from == restart_water
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=3).status == :packing_ready
    @test sys3.molecule_positions[1].cm ≈ saved_water_cm
    rm(sys3.output_file; force=true)

    rm(restart_all; force=true)
    rm(restart_water; force=true)
end

@testitem "restart_from a PDB file" begin
    using PDBTools: read_pdb

    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    # First run: a normal PDB output file, no restart involved — it's not
    # Packmol's own restart format, but it should work as a restart_from
    # source anyway, since it's just atom positions like any other PDB.
    water1 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys1 = PackmolSystem([water1]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3))
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=11).status == :packing_ready
    source_pdb = sys1.output_file

    # Restarting from that PDB (file path, dispatched by its .pdb extension)
    # should reproduce essentially the same, already-converged positions —
    # converging immediately — using a different seed to rule out the run
    # just coincidentally reconverging to the same configuration on its own.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3),
        restart_from=source_pdb,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=22).status == :packing_ready
    for i in 1:4
        @test sys2.molecule_positions[i].cm ≈ sys1.molecule_positions[i].cm atol = 1e-2
    end
    rm(sys2.output_file; force=true)

    # Same, but passing an already-loaded `Vector{<:Atom}` directly (the
    # Julia-API-only form) instead of a file path.
    atoms = read_pdb(source_pdb)
    water3 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([water3]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3),
        restart_from=atoms,
    )
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=33).status == :packing_ready
    for i in 1:4
        @test sys3.molecule_positions[i].cm ≈ sys1.molecule_positions[i].cm atol = 1e-2
    end
    rm(sys3.output_file; force=true)

    rm(source_pdb; force=true)
end

@testitem "restart_from with a fixed molecule present" begin
    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    restart_all = tempname()

    # Fixed structure type declared *before* the free one — matching common
    # real-world inputs (e.g. a fixed protein listed first, solvent after) —
    # so a whole-system restart must correctly skip over the fixed
    # molecule's own atoms/lines (interleaved at the front of the restart
    # source) instead of assuming the source only ever lists free molecules.
    fixed_water = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water1 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    # `avoid_overlap=false`: with only 3 fixed atoms, the fixed-particle
    # overlap-check system's own auto-sized cutoff (see
    # `_build_overlap_check_system`) can exceed half its own bounding box —
    # a pre-existing limitation for tiny fixed structures, unrelated to what
    # this test is actually checking (restart_from's atom bookkeeping).
    sys1 = PackmolSystem([fixed_water, water1]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_to=restart_all, avoid_overlap=false,
    )
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=41).status == :packing_ready
    source_pdb = sys1.output_file
    saved_free_cm = [sys1.molecule_positions[i].cm for i in 2:5]

    # Whole-system restart_from a PDB written for the entire system (fixed
    # molecule included) must not error on the atom count, and must recover
    # the free molecules' positions (the fixed one keeps its own `fixed`
    # position regardless of what the restart source contains for it).
    fixed_water2 = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([fixed_water2, water2]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=source_pdb, avoid_overlap=false,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=42).status == :packing_ready
    for i in 2:5
        @test sys2.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys2.output_file; force=true)

    # Same, but restarting from the raw whole-system restart_to file (also
    # written for the entire system, fixed molecule included).
    fixed_water3 = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water3 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([fixed_water3, water3]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all, avoid_overlap=false,
    )
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=43).status == :packing_ready
    for i in 2:5
        @test sys3.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys3.output_file; force=true)

    rm(source_pdb; force=true)
    rm(restart_all; force=true)
end
