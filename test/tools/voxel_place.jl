@testitem "voxel_place" setup=[TestFiles] begin
    using Packmol
    using StaticArrays: SVector
    water = TestFiles.test_file("structure_files/water.pdb")
    protein = TestFiles.test_file("structure_files/protein.pdb")
    quiet(f) = redirect_stdout(f, devnull)
    output = tempname() * ".pdb"

    # Everything fits: all placed, starting from their positions, and packed from them
    sts = [structure_type(protein; number = 3, constraints = [InsideSphere([0., 0., 0.], 80.)]),
           structure_type(water; number = 1000, constraints = [InsideSphere([0., 0., 0.], 80.)])]
    sys = PackmolSystem(sts; output, tolerance = 2.0)
    placed, dropped = quiet(() -> voxel_place(sys))
    @test dropped == [0, 0]
    @test [st.number_of_molecules for st in placed.structure_types] == [3, 1000]
    @test placed.restart_from isa Vector
    @test length(placed.restart_from) == length(placed.atoms)
    @test sys.status == :not_packed
    # the placed atoms are within the constraints, and farther than the tolerance
    x = [SVector{3}(a.x, a.y, a.z) for a in placed.restart_from]
    @test all(xi -> sqrt(sum(abs2, xi)) < 80.0 + 0.1, x)
    mol = [a.molecule_index for a in placed.atoms]
    dmin = minimum(sqrt(sum(abs2, x[i] - x[j])) for i in eachindex(x), j in eachindex(x) if mol[i] < mol[j])
    @test dmin > 2.0
    @test quiet(() -> packmol(placed)).status == :packing_ready

    # Not everything fits: the molecules that do not fit are dropped, and the rest pack
    sts = [structure_type(protein; number = 30, constraints = [InsideSphere([0., 0., 0.], 70.)])]
    sys = PackmolSystem(sts; output, tolerance = 2.0)
    result = Ref{Any}()
    log = mktemp() do path, io
        redirect_stdout(io) do
            result[] = voxel_place(sys; max_tries_mult = 100)
        end
        close(io)
        read(path, String)
    end
    placed, dropped = result[]
    @test occursin("SIGNIFICANT WARNING", log)
    n = placed.structure_types[1].number_of_molecules
    @test 0 < n && dropped[1] > 0
    @test n + dropped[1] == 30
    @test quiet(() -> packmol(placed)).status == :packing_ready

    # Periodic cell, with a fixed molecule, from an input file
    dir = mktempdir()
    cp(water, joinpath(dir, "water.pdb"))
    cp(protein, joinpath(dir, "protein.pdb"))
    write(joinpath(dir, "box.inp"), """
        tolerance 2.0
        pbc 120. 120. 120.
        output box.pdb
        structure protein.pdb
          number 1
          fixed 60. 60. 60. 0. 0. 0.
        end structure
        structure water.pdb
          number 4000
        end structure
        """)
    placed, dropped = quiet(() -> voxel_place(joinpath(dir, "box.inp")))
    @test dropped == [0, 0]
    @test quiet(() -> packmol(placed)).status == :packing_ready

    # The input keyword: the voxel placement, then the packing of the molecules that fit
    write(joinpath(dir, "box_voxel.inp"), "voxel_place 500\n" * read(joinpath(dir, "box.inp"), String))
    sys = Packmol.read_packmol_input(joinpath(dir, "box_voxel.inp"))
    @test sys.voxel_place == 500
    sys = quiet(() -> packmol(sys))
    @test sys.status == :packing_ready
    @test sys.restart_from isa Vector
    @test Packmol.read_packmol_input(joinpath(dir, "box.inp")).voxel_place === nothing
    write(joinpath(dir, "box_voxel.inp"), "voxel_place\n" * read(joinpath(dir, "box.inp"), String))
    @test Packmol.read_packmol_input(joinpath(dir, "box_voxel.inp")).voxel_place == 1000

    # 3D only, and not from a restart
    sys = PackmolSystem([structure_type(water; number = 2)]; output, restart_from = output)
    @test_throws ArgumentError voxel_place(sys)
    rm(output; force = true)
end
