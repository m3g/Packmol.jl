@testitem "structure_type" setup=[TestFiles] begin
    using StaticArrays
    file = TestFiles.test_file("structure_files/water.pdb")
    tolerance = 2.0

    st = structure_type(file; number=1000, constraints=[InsideBox([0.,0.,0.],[40.,40.,40.])])
    @test st.filename == file
    @test st.natoms == 3
    @test st.number_of_molecules == 1000
    @test st.fixed.fixed == false
    @test st.constraints == [Box{Inside,Float64}([20.,20.,20.],[40.,40.,40.],5.0)]
    @test st.atom_constraints == [[1], [1], [1]]
    @test st.radii == [1.0, 1.0, 1.0]  # defaults to tolerance/2

    st2 = structure_type(file; number=1, radius=0.5)
    @test st2.radii == [0.5, 0.5, 0.5]

    # fixed + geometric center
    st3 = structure_type(file; number=1, fixed=([1.0,2.0,3.0],[0.0,0.0,0.0]), center=:geometric)
    @test st3.fixed.fixed == true
    @test st3.fixed.position.cm == SVector(1.0,2.0,3.0)
    @test sum(st3.reference_coordinates) / length(st3.reference_coordinates) ≈ SVector(1.0,2.0,3.0)

    # center without fixed is an error
    @test_throws ArgumentError structure_type(file; number=1, center=:geometric)

    # T=Float32
    st4 = structure_type(file; number=1, T=Float32)
    @test eltype(st4.radii) == Float32
    @test eltype(st4.reference_coordinates) == SVector{3,Float32}
end

#=
    read_structure_data(input_file_block::IOBuffer, tolerance; T=Float64, D=3)
    read_structure_data(input_file_block::AbstractString, args... kargs..)

Reads the structure data from a structure block of the input file. Requires the `tolerance`
parameter to set atom radii by default. The function returns a `StructureType` object. 

The type `T` defines the number type (Float32, Float64, etc.) and `D` the number of dimensions (2 or 3).

=#

@testitem "read_structure_data" setup=[RigidBody, TestFiles] begin
    using Packmol
    using PDBToolsBase: read_pdb, positions

    file = TestFiles.test_file("structure_files/water.pdb")
    tolerance = 2.0

    input_file_block = """
    structure $file        
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test isfile(s.filename)
    @test length(s.atoms) == 3
    @test s.number_of_molecules == 1000
    @test length(s.reference_coordinates) == 3
    @test length(s.constraints) == 1
    # "inside box" stores (center, sides) computed from the given (lo, hi) corners:
    # a box from (0,0,0) to (40,40,40) has center (20,20,20).
    @test s.constraints[1] == Box{Inside, Float64}([20.0, 20.0, 20.0], [40.0, 40.0, 40.0], 5.0)
    @test s.fixed.fixed == false
    @test s.atom_constraints == [[1], [1], [1]]
    # Default per-atom radius is tolerance/2, so that two touching atoms'
    # radii sum to the requested tolerance (see cartesian_fg!).
    @test s.radii == [1.0, 1.0, 1.0]

    input_file_block = """
    structure $file
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
        outside sphere 0. 0. 0. 10.
        radius 1.0
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test s.radii == [1.0, 1.0, 1.0]
    @test s.atom_constraints == [[1, 2], [1, 2], [1, 2]]
    @test s.constraints[1] == Box{Inside, Float64}([20.0, 20.0, 20.0], [40.0, 40.0, 40.0], 5.0)
    # Sphere's default constraint weight (1e-2) differs from Box's (5.0).
    @test s.constraints[2] == Sphere{Outside, Float64}([0.0, 0.0, 0.0], 10.0, 1e-2)

    input_file_block = """
    structure $file
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
        atoms 1 3
            outside sphere 0. 0. 0. 10.
        end
        atoms 1 2
            radius 1.0
        end
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test s.radii == [1.0, 1.0, 1.0]
    @test s.atom_constraints == [[1, 2], [1], [1, 2]]
    @test s.constraints[1] == Box{Inside, Float64}([20.0, 20.0, 20.0], [40.0, 40.0, 40.0], 5.0)
    @test s.constraints[2] == Sphere{Outside, Float64}([0.0, 0.0, 0.0], 10.0, 1e-2)

    s = Packmol.read_structure_data(input_file_block, tolerance; T = Float32)
    @test s.radii == Float32[1.0, 1.0, 1.0]
    @test s.constraints[1] == Box{Inside, Float32}([20.0, 20.0, 20.0], [40.0, 40.0, 40.0], 5.0)
    @test s.constraints[2] == Sphere{Outside, Float32}([0.0, 0.0, 0.0], 10.0, 1e-2)


end

@testitem "fixed molecules" setup=[RigidBody, TestFiles] begin
    using Packmol
    using PDBToolsBase: positions
    file = TestFiles.test_file("structure_files/diatomic.pdb")
    tolerance = 2.0

    # Fixed molecule: do not move
    input_file_block = """
    structure $file        
        number 1
        fixed 0. 0. 0. 0. 0. 0.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test s.fixed.fixed == true
    @test positions(s.atoms) ≈ s.reference_coordinates

    # Fixed molecule without rotation: center of mass at origin
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 0. 0. 0. 0. 0. 0.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]] 

    # Fixed molecule without rotation: center of mass at (10, 10, 10)
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 10. 10. 10. 0. 0. 0.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[10.0, 10.0, 9.5], [10.0, 10.0, 10.5]]

    # Rotate fixed molecule π/2 around x-axis (counterclockwise)
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 0. 0. 0. $(π/2) 0. 0.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]]

    # Rotate fixed molecule π/2 around y-axis (counterclockwise)
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 0. 0. 0. 0. $(π/2) 0. 
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]

    # Rotate fixed molecule π/2 around z-axis (counterclockwise)
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 0. 0. 0. 0. 0. $(π/2)
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]]

    # Rotate fixed molecule π/2 around z-axis (counterclockwise) and move it to (10, 10, 10)
    input_file_block = """
    structure $file        
        number 1
        center
        fixed 10. 10. 10. 0. 0. $(π/2)
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance) 
    @test s.reference_coordinates ≈ SVector{3, Float64}[[10.0, 10.0, 9.5], [10.0, 10.0, 10.5]]

end

@testitem "structure_type: atom_constraints" setup=[TestFiles] begin
    file = TestFiles.test_file("structure_files/palmitoil.pdb")
    st = structure_type(file;
        number = 2,
        constraints = [InsideBox([0., 0., 0.], [40., 40., 14.])],
        atom_constraints = [
            [31, 32] => BelowPlane([0., 0., 1.], 2.0),
            [1, 2] => [AbovePlane([0., 0., 1.], 12.0)],
        ],
    )
    # The same as the `atoms ... end atoms` blocks of an input file
    st_input = Packmol.read_structure_data("""
        structure $file
          number 2
          inside box 0. 0. 0. 40. 40. 14.
          atoms 31 32
            below plane 0. 0. 1. 2.
          end atoms
          atoms 1 2
            above plane 0. 0. 1. 12.
          end atoms
        end structure
        """, 2.0)
    @test st.constraints == st_input.constraints
    @test st.atom_constraints == st_input.atom_constraints
    @test st.atom_constraints[31] == [1, 2]
    @test st.atom_constraints[1] == [1, 3]
    @test st.atom_constraints[10] == [1]
    @test_throws ArgumentError structure_type(file; number=1, atom_constraints=[[1000] => BelowPlane([0., 0., 1.], 2.0)])
end

@testitem "structure_type: several conformations" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb
    water = TestFiles.test_file("structure_files/water.pdb")
    lines = filter(l -> startswith(l, "ATOM") || startswith(l, "HETATM"), readlines(water))
    dir = mktempdir()
    # Three models of water (the same coordinates: only the models matter here)
    models = joinpath(dir, "water_models.pdb")
    open(models, "w") do io
        for k in 1:3
            println(io, "MODEL     ", lpad(k, 4))
            foreach(l -> println(io, l), lines)
            println(io, "ENDMDL")
        end
        println(io, "END")
    end
    sts = structure_type(models; number=7, constraints=[InsideBox([0., 0., 0.], [20., 20., 20.])])
    @test sts isa Vector{<:Packmol.StructureType}
    @test [st.number_of_molecules for st in sts] == [3, 2, 2]
    @test [st.conformation for st in sts] == [(1, 3), (2, 3), (3, 3)]
    @test all(st.natoms == length(lines) for st in sts)
    @test all(st.constraints == sts[1].constraints && st.constraints !== sts[1].constraints for st in sts[2:end])
    # Fewer molecules than conformations
    sts2 = structure_type(models; number=2)
    @test [st.conformation for st in sts2] == [(1, 2), (2, 2)]
    # A vector of files
    sts3 = structure_type([water, water]; number=4)
    @test [st.number_of_molecules for st in sts3] == [2, 2]
    # A single conformation is a single structure type
    @test structure_type(water; number=4) isa Packmol.StructureType
    # Errors: different atoms, fixed molecules, restart files
    other = TestFiles.test_file("structure_files/palmitoil.pdb")
    @test_throws ArgumentError structure_type([water, other]; number=2)
    @test_throws ArgumentError structure_type(models; number=1, fixed=([0., 0., 0.], [0., 0., 0.]))
    @test_throws ArgumentError structure_type(models; number=3, restart_to=joinpath(dir, "r.pdb"))

    # The same in an input file, packed: the conformations are a single structure type
    # in the show of the system, and the molecules are numbered as a single type
    write(joinpath(dir, "conformations.inp"), """
        tolerance 2.0
        output conformations.pdb
        structure water_models.pdb
          number 7
          inside box 0. 0. 0. 20. 20. 20.
        end structure
        structure $water
          number 5
          inside box 0. 0. 0. 20. 20. 20.
        end structure
        """)
    sys = Packmol.read_packmol_input(joinpath(dir, "conformations.inp"))
    @test length(sys.structure_types) == 4
    @test sys.nmols == 12
    shown = sprint(show, MIME"text/plain"(), sys)
    @test occursin("Number of types of structures: 2", shown)
    @test occursin("\"water_models.pdb\" - 3 conformations", shown)
    @test occursin("n_mols = 7", shown)
    redirect_stdout(devnull) do
        packmol(sys)
    end
    @test sys.status == :packing_ready
    out = read_pdb(joinpath(dir, "conformations.pdb"))
    resnums = [at.resnum for at in out if at.name == "OH2" || at.name == "O"]
    @test resnums[1:7] == 1:7
    @test resnums[8:12] == 1:5
    # PackmolSystem from a vector of structure types and vectors of them
    sys2 = PackmolSystem([sts, structure_type(water; number=5)]; output="")
    @test length(sys2.structure_types) == 4
end

@testitem "shape_descriptors and select_conformations" setup=[TestFiles] begin
    using StaticArrays: SVector
    # A rod: atoms on a line (anisotropy 1); a sphere-like set: the vertices of a cube
    rod = [SVector(0.0, 0.0, z) for z in 0.0:1.0:10.0]
    d = shape_descriptors(rod)
    @test d.anisotropy ≈ 1.0
    @test d.asphericity ≈ d.principal_variances[1]
    @test d.acylindricity ≈ 0.0 atol = 1e-12
    cube = [SVector(x, y, z) for x in (-1.0, 1.0), y in (-1.0, 1.0), z in (-1.0, 1.0)][:]
    d = shape_descriptors(cube)
    @test d.anisotropy ≈ 0.0 atol = 1e-12
    @test d.radius_of_gyration ≈ sqrt(3.0)
    disc = [SVector(cos(θ), sin(θ), 0.0) for θ in range(0, 2π; length=13)[1:12]]
    @test shape_descriptors(disc).anisotropy ≈ 0.25

    # A file with three models of water, stretched along z by different amounts
    water = TestFiles.test_file("structure_files/water.pdb")
    lines = filter(l -> startswith(l, "HETATM") || startswith(l, "ATOM"), readlines(water))
    dir = mktempdir()
    models = joinpath(dir, "models.pdb")
    open(models, "w") do io
        for (k, s) in enumerate((1.0, 3.0, 2.0))
            println(io, "MODEL     ", lpad(k, 4))
            for l in lines
                z = parse(Float64, l[47:54])
                y = parse(Float64, l[39:46])
                println(io, l[1:38], lpad(string(round(s * y; digits=3)), 8), lpad(string(round(s * z; digits=3)), 8), l[55:end])
            end
            println(io, "ENDMDL")
        end
    end
    ds = shape_descriptors(models)
    @test length(ds) == 3
    @test shape_descriptors(water) isa NamedTuple
    out = joinpath(dir, "rods.pdb")
    indices, descriptors = select_conformations(models, out; n=2)
    @test indices == [2, 3]
    @test [d.asphericity for d in descriptors] == [ds[2].asphericity, ds[3].asphericity]
    sts = structure_type(out; number=4)
    @test length(sts) == 2
    @test shape_descriptors(out)[1].asphericity ≈ ds[2].asphericity rtol = 1e-3
    indices, _ = select_conformations(models, joinpath(dir, "rod.pdb")) # n = 1 by default
    @test indices == [2]
    @test structure_type(joinpath(dir, "rod.pdb"); number=4) isa Packmol.StructureType
    indices, _ = select_conformations(models, joinpath(dir, "spheres.pdb"); n=1, shape=:sphere)
    @test length(indices) == 1
    @test_throws ArgumentError select_conformations(models, out; n=1, shape=:disc)
end
