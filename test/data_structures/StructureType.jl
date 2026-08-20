@testitem "structure_type" begin
    using StaticArrays
    file = Packmol.src_dir * "/../test/structure_files/water.pdb"
    tolerance = 2.0

    st = structure_type(file; number=1000, constraints=[InsideBox([0.,0.,0.],[40.,40.,40.])])
    @test st.filename == file
    @test st.natoms == 3
    @test st.number_of_molecules == 1000
    @test st.fixed.fixed == false
    @test st.constraints == [Box{Inside,Float64}([0.,0.,0.],[40.,40.,40.],5.0)]
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

@testitem "read_structure_data" setup=[RigidBody] begin
    using Packmol
    using PDBTools

    file = Packmol.src_dir*"/../test/structure_files/water.pdb"
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

@testitem "fixed molecules" setup=[RigidBody] begin
    using Packmol, PDBTools
    file = Packmol.src_dir*"/../test/structure_files/diatomic.pdb"
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
    @test coor(s.atoms) ≈ s.reference_coordinates

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
