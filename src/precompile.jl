#
# Precompilation workload: tiny systems (a few molecules, one loop) that go through the
# input parser, the constraints, periodic cells, the packing loop and the output formats,
# so that the first packing of a Julia session does not wait for their compilation.
#
# Not run when building the executable (juliac/build.jl sets
# PACKMOL_PRECOMPILE_WORKLOAD=false): `juliac --trim` would reuse the code compiled
# here, and the trim verifier rejects some of it (error paths of SPGBox that build
# their messages dynamically), while the same functions compiled for the executable pass.
#
using PrecompileTools: @setup_workload, @compile_workload

if get(ENV, "PACKMOL_PRECOMPILE_WORKLOAD", "true") != "false"
@setup_workload begin
    data = joinpath(RecipesDirectory, "test", "data")
    water = joinpath(data, "water.pdb")
    ethanol = joinpath(data, "ethanol.pdb")
    popc = joinpath(RecipesDirectory, "lipids", "popc.pdb")
    @compile_workload begin
        dir = mktempdir()
        redirect_stdio(stdout=devnull, stderr=devnull) do
            # An input file, without periodic boundary conditions, with the constraints
            # of the input syntax, a fixed molecule, constraints on some atoms and a
            # rotation constraint
            input = joinpath(dir, "constraints.inp")
            write(input, """
            tolerance 2.0
            output constraints.pdb
            nloop 1
            seed 1
            structure $ethanol
              number 1
              fixed 0. 0. 0. 0. 0. 0.
              center
            end structure
            structure $water
              number 2
              inside box -20. -20. -20. 20. 20. 20.
              outside sphere 0. 0. 0. 4.
            end structure
            structure $water
              number 2
              inside cube -20. -20. -20. 40.
              outside box -3. -3. -3. 3. 3. 3.
              outside cube -3. -3. -3. 6.
            end structure
            structure $water
              number 2
              inside sphere 0. 0. 0. 20.
              above plane 0. 0. 1. -15.
              below plane 0. 0. 1. 15.
            end structure
            structure $water
              number 2
              inside cylinder 0. 0. -15. 0. 0. 1. 15. 30.
              outside cylinder 0. 0. -15. 0. 0. 1. 3. 30.
            end structure
            structure $water
              number 2
              inside ellipsoid 0. 0. 0. 20. 15. 15. 1.
              outside ellipsoid 0. 0. 0. 4. 4. 4. 1.
            end structure
            structure $water
              number 2
              inside box -20. -20. -20. 20. 20. 20.
              above gaussian 0. 0. 1. 1. 0. 0. -15. 2. 0. 3.
              below sin 0. 0. 1. 1. 0. 0. 15. 2. 10. 0.
              above exponential 0. 1. 0. 1. 0. 0. -15. 1. 0. 0.1
              below radial_gaussian 0. 0. 1. 0. 0. 0. 15. 2. 3.
            end structure
            structure $ethanol
              number 2
              resnumbers 2
              inside box -20. -20. -20. 20. 20. 20.
              atoms 1
                radius 1.5
                above plane 0. 0. 1. -10.
              end atoms
              constrain_rotation z 0. 30.
                atoms 1
                to 6
              end constrain_rotation
            end structure
            structure $popc
              number 1
              inside box -20. -20. -20. 20. 20. 20.
              constrain_rotation z 0. 20.
            end structure
            """)
            packmol(input; iprint = 1000)
            # Checking the initial placement, and an input with an orthorhombic periodic
            # box and XYZ output
            sys = read_packmol_input(input)
            sys.check_initial_placement = true
            packmol(sys)
            input_pbc = joinpath(dir, "pbc.inp")
            write(input_pbc, """
            tolerance 2.0
            output pbc.xyz
            nloop 1
            pbc 15. 15. 15.
            structure $water
              number 5
            end structure
            structure $ethanol
              number 2
              above plane 0. 0. 1. 5.
            end structure
            """)
            packmol(input_pbc; iprint = 1000)
            free_volume(input_pbc; npoints = 1000)
            guess_volume(input_pbc)
            # The Julia API: a rhombic dodecahedron, written in its shape, as mmCIF
            w = structure_type(water; number = 5, constraints = [OutsideSphere([0., 0., 0.], 3.)])
            e = structure_type(ethanol; number = 2,
                atom_constraints = [[1] => AbovePlane([0., 0., 1.], 0.)],
                constrain_rotation = Dict(:z => (0.0, 30.0, [1], [6])),
            )
            sys = PackmolSystem([w, e]; tolerance = 2.0, output = joinpath(dir, "julia.cif"),
                unitcell = dodecahedral_unitcell(Float64, 16.0), nloop = 1)
            triclinic_to_dodecahedral(sys)
            packmol(sys; iprint = 1000)
            get_atoms(sys)
            sys = PackmolSystem([w]; tolerance = 2.0, output = joinpath(dir, "octahedral.pdb"),
                unitcell = octahedral_unitcell(Float64, 16.0), nloop = 1)
            packmol(sys; iprint = 1000)
            # A unit cell given as a matrix, and a rotation constraint of the longest axis
            e = structure_type(ethanol; number = 2, constraints = [InsideBox([0., 0., 0.], [15., 15., 15.])],
                constrain_rotation = Dict(:z => (0.0, 30.0)))
            sys = PackmolSystem([w, e]; tolerance = 2.0, output = joinpath(dir, "matrix.pdb"),
                unitcell = [15. 0. 0.; 0. 15. 0.; 0. 0. 15.], nloop = 1)
            packmol(sys; iprint = 1000)
            # A recipe (the system is not packed: it has thousands of molecules)
            recipe = SolutionBoxUS(solute_pdbfile = joinpath(data, "poly_h.pdb"),
                solvent_pdbfile = water, density = 1.0u"g/mL")
            show(stdout, MIME"text/plain"(), recipe)
            write_packmol_input(recipe; margin = 20.0, input = joinpath(dir, "recipe.inp"),
                output = joinpath(dir, "recipe.pdb"))
        end
    end
end
end
