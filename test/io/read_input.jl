@testitem "_parse_value" begin
    using Packmol: _parse_value
    @test _parse_value(Int, "max_iter", "100") == 100
    @test _parse_value(Float64, "tolerance", "2.0") == 2.0
    @test _parse_value(Float32, "tolerance", "2.0") == 2.0f0
    @test_throws ArgumentError _parse_value(Int, "max_iter", "100.0")
end

@testitem "_parse_options" begin
    using Packmol: _parse_options
    @test _parse_options("connect", "yes", ("yes" => true, "no" => false)) == true
    @test _parse_options("connect", "no", ("yes" => true, "no" => false)) == false
    @test_throws ArgumentError _parse_options("connect", "nop", ("yes" => true, "no" => false))
end

@testitem "read_packmol_input" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    file = TestFiles.test_file("input_files/water_box.inp")
    sys = read_packmol_input(file)
    @test sys.filetype == "pdb"
    @test sys.output_file == joinpath(dirname(file), "water_box.pdb")
    @test sys.nmols == 1000
    @test length(sys.atoms) == 3000
    @test all(at.molecule_index == 1 for at in sys.atoms[1:3])
    @test all(at.molecule_index == 2 for at in sys.atoms[4:6])
    @test all(at.molecule_index == 1000 for at in sys.atoms[2998:3000])
    @test all(at.radius ≈ 1.0 for at in sys.atoms)
    @test length(sys.molecule_positions) == 1000
    @test isnothing(sys.radscale)
    # default scales: 1.0 for water (3 atoms)
    @test Packmol._radscale(sys, sys.structure_types[1]) == 1.0
    # fixed molecules: the largest scale of the free ones (5.0 with a lipid, 1.0 with water only)
    water = TestFiles.test_file("structure_files/water.pdb")
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    fx = structure_type(water; number=1, fixed=([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    w = structure_type(water; number=2, constraints=[InsideBox([0, 0, 0], [20, 20, 20])])
    l = structure_type(lipid; number=2, constraints=[InsideBox([0, 0, 0], [20, 20, 20])])
    s2 = PackmolSystem([fx, w, l]; output=tempname() * ".pdb")
    @test [Packmol._radscale(s2, st) for st in s2.structure_types] == [5.0, 1.0, 5.0]
    s2 = PackmolSystem([fx, w]; output=tempname() * ".pdb")
    @test [Packmol._radscale(s2, st) for st in s2.structure_types] == [1.0, 1.0]
    s2 = PackmolSystem([fx, w, structure_type(lipid; number=2, radscale=2.0)]; output=tempname() * ".pdb")
    @test Packmol._radscale(s2, s2.structure_types[1]) == 2.0
    # restarted molecules: 1.0 (the system or the type restarted), unless set
    s2 = PackmolSystem([fx, w, l]; output=tempname() * ".pdb")
    s2.restart_from = "restart.pdb"
    @test [Packmol._radscale(s2, st) for st in s2.structure_types] == [1.0, 1.0, 1.0]
    lr = structure_type(lipid; number=2, restart_from="restart.pdb", constraints=[InsideBox([0, 0, 0], [20, 20, 20])])
    s2 = PackmolSystem([fx, w, lr, l]; output=tempname() * ".pdb")
    @test [Packmol._radscale(s2, st) for st in s2.structure_types] == [5.0, 1.0, 1.0, 5.0]
    s2 = PackmolSystem([fx, w, structure_type(lipid; number=2, restart_from="restart.pdb", radscale=3.0)]; output=tempname() * ".pdb")
    @test [Packmol._radscale(s2, st) for st in s2.structure_types] == [3.0, 1.0, 3.0]

    sys = read_packmol_input(file; T=Float32)
    @test typeof(sys.tolerance) == Float32
    @test eltype(sys.molecule_positions) == Packmol.MoleculePosition{3,Float32}
    @test isnothing(sys.radscale)

    # 2D: Currently not supported
    # sys = read_packmol_input(file; D=2)
end

@testitem "radscale / discale keyword" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    dir = dirname(TestFiles.test_file("input_files/water_box.inp"))
    original = read(TestFiles.test_file("input_files/water_box.inp"), String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))

        radscale_file = joinpath(tmp, "radscale.inp")
        write(radscale_file, original * "\nradscale 1.5\n")
        sys = read_packmol_input(radscale_file)
        @test sys.radscale == 1.5
        @test isnothing(sys.radscale_decay_loops)

        decay_file = joinpath(tmp, "radscale_decay.inp")
        write(decay_file, original * "\nradscale 3.0\nradscale_decay_loops 50\n")
        sys = read_packmol_input(decay_file)
        @test sys.radscale == 3.0
        @test sys.radscale_decay_loops == 50

        discale_file = joinpath(tmp, "discale.inp")
        write(discale_file, original * "\ndiscale 1.3\n")
        sys = @test_warn r"discale" read_packmol_input(discale_file)
        @test sys.radscale == 1.3
    end
end

@testitem "check_initial_placement / check keyword" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    dir = dirname(TestFiles.test_file("input_files/water_box.inp"))
    original = read(TestFiles.test_file("input_files/water_box.inp"), String)
    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))
        file = joinpath(tmp, "default.inp")
        write(file, original)
        @test !read_packmol_input(file).check_initial_placement
        write(file, original * "\ncheck_initial_placement\n")
        @test read_packmol_input(file).check_initial_placement
        # Commented, as in the tutorials
        write(file, original * "\n# check_initial_placement\n")
        @test !read_packmol_input(file).check_initial_placement
        # Fortran's name
        write(file, original * "\ncheck\n")
        sys = @test_warn r"check_initial_placement" read_packmol_input(file)
        @test sys.check_initial_placement
    end
end

@testitem "per-structure keyword does not leak into the global setting" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    # `restart_from` is valid both as a global keyword and as a
    # per-structure one (inside `structure ... end structure`); a
    # per-structure occurrence must only set that structure type's field,
    # not the whole-system `PackmolSystem` field of the same name —
    # regression test for a parser bug where the global-keyword dispatch
    # table was consulted unconditionally, before checking whether the
    # current line was inside a structure block.
    dir = dirname(TestFiles.test_file("input_files/water_box.inp"))
    original = read(TestFiles.test_file("input_files/water_box.inp"), String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))
        restart_touched = joinpath(tmp, "some_restart_file.pdb")

        text = replace(original, "end structure" => "restart_from $restart_touched\nend structure")
        file = joinpath(tmp, "per_structure_restart.inp")
        write(file, text)

        sys = read_packmol_input(file)
        @test isnothing(sys.restart_from)
        @test only(sys.structure_types).restart_from == restart_touched
    end
end

@testitem "read_packmol_input: comments, legacy keywords, nloop" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    input = joinpath(dir, "test.inp")
    write(input, """
    # A comment line
    tolerance 2.0 # inline comment
    output out.pdb
    nloop 37
    maxit 15
    movefrac 0.1
    nthreads 2
    movebadrandom
    short_tol_dist 1.0
    sidemax 1000.
    writecrd out.crd
    randominitialpoint
    optimization_print_level 2
    chkgrad
    writeout 10
    writebad
    structure $water # the water
      number 10
      nloop 5
      segid WAT
      inside box 0. 0. 0. 20. 20. 20. # inline comment
    end structure
    """)
    warnings = joinpath(dir, "warnings.txt")
    sys = open(warnings, "w") do io
        redirect_stderr(io) do
            read_packmol_input(input)
        end
    end
    for keyword in ("randominitialpoint", "optimization_print_level", "chkgrad", "writebad")
        @test occursin("$keyword is a deprecated keyword and was ignored", read(warnings, String))
    end
    @test sys.tolerance == 2.0
    @test sys.nloop == 37
    @test sys.max_iter == 15
    @test sys.movefrac == 0.1
    @test sys.nthreads == 2
    @test sys.writeout == 10
    @test sys.structure_types[1].segid == "WAT"
    @test sys.nmols == 10
    @test length(sys.structure_types[1].constraints) == 1
    # Not set in the input file: the `packmol` defaults apply
    sys = read_packmol_input(TestFiles.test_file("input_files/water_box.inp"))
    @test isnothing(sys.nloop)
    @test isnothing(sys.max_iter)
    @test isnothing(sys.movefrac)
    @test isnothing(sys.nthreads)
    @test isnothing(sys.structure_types[1].segid)
    # segid has at most 4 characters
    write(input, "tolerance 2.0\nstructure $water\n  number 1\n  segid TOOLONG\nend structure\n")
    @test_throws ArgumentError read_packmol_input(input)
    # Misspelled keywords in structure blocks are errors; deprecated and global keywords
    # are ignored there (with a warning)
    write(input, "tolerance 2.0\nstructure $water\n  numbr 10\nend structure\n")
    @test_throws ArgumentError read_packmol_input(input)
    write(input, "tolerance 2.0\nstructure $water\n  number 10\n  seed 1\n  fscale 1.0\n  inside box 0. 0. 0. 20. 20. 20.\nend structure\n")
    sys = redirect_stderr(devnull) do
        read_packmol_input(input)
    end
    @test sys.nmols == 10
end

@testitem "read_packmol_input from IO" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    # Paths are relative to the current directory
    dir = TestFiles.test_file("structure_files")
    input = "tolerance 2.0\noutput out.pdb\nstructure water.pdb\n  number 10\n  inside box 0. 0. 0. 20. 20. 20.\nend structure\n"
    sys = cd(dir) do
        read_packmol_input(IOBuffer(input))
    end
    @test sys.input_file == ""
    @test sys.output_file == "out.pdb"
    @test sys.nmols == 10
    @test sys.structure_types[1].filename == "water.pdb"
end

@testitem "read_packmol_input: constrain_rotation" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    file = tempname() * ".inp"
    write(file, """
        tolerance 2.0
        output $(tempname()).pdb
        structure $lipid
          number 10
          inside box -20 -20 -20 20 20 20
          constrain_rotation z 0. 15.
            atoms 1 2   # head
            to 42
          end constrain_rotation
          constrain_rotation x 90. 10.
          atoms 42
            below plane 0. 0. 1. 15.
          end atoms
        end structure
        """)
    sys = read_packmol_input(file)
    st = only(sys.structure_types)
    rz, rx = st.rotation_constraints
    @test rz.axis == [0, 0, 1] && rz.from == [1, 2] && rz.to == [42]
    @test rz.lower == 0 && rz.upper ≈ deg2rad(15)
    @test rx.axis == [1, 0, 0] && isempty(rx.from) && rx.lower ≈ deg2rad(80) && rx.upper ≈ deg2rad(100)
    # The `atoms` block of per-atom constraints is still read
    @test length(st.constraints) == 2
    @test st.atom_constraints[42] == [1, 2] && st.atom_constraints[1] == [1]
    # Errors
    for (text, msg) in (
        ("constrain_rotation w 0. 15.", "axis"),
        ("constrain_rotation z 0. 15.\n atoms 1\n to 99\n end constrain_rotation", "out of range"),
        ("constrain_rotation z 0. 15.\n atoms 1\n to 2", "end constrain_rotation"),
        ("fixed 0. 0. 0. 0. 0. 0.\n constrain_rotation z 0. 15.", "fixed"),
    )
        write(file, "tolerance 2.0\noutput $(tempname()).pdb\nstructure $lipid\n number 1\n $text\nend structure\n")
        err = try
            read_packmol_input(file); nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin(msg, err.msg)
    end
end

@testitem "read_packmol_input: periodic_boundary_style" setup=[TestFiles] begin
    water = TestFiles.test_file("structure_files/water.pdb")
    input(lines) = IOBuffer("tolerance 2.0\noutput out.pdb\n$lines\nstructure $water\n  number 10\nend structure\n")
    @test read_packmol_input(input("pbc dodecahedral 0. 0. 0. 20.")).periodic_boundary_style == :triclinic
    for style in (:triclinic, :dodecahedral, :octahedral)
        sys = read_packmol_input(input("pbc dodecahedral 0. 0. 0. 20.\nperiodic_boundary_style $style"))
        @test sys.periodic_boundary_style == style
    end
    @test_throws ArgumentError read_packmol_input(input("pbc 20. 20. 20.\nperiodic_boundary_style compact"))
    @test_throws ArgumentError read_packmol_input(input("periodic_boundary_style dodecahedral"))
end
