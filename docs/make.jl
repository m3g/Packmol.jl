using Documenter, DocumenterCodeBlocks
using Markdown
using Packmol
push!(LOAD_PATH, "../src/")

# The input files of the tutorials are shown in the pages from the files themselves
# (docs/src/tutorials/files), with an `@eval` block (in a drop-down box, between
# `<details>` tags of `@raw html` blocks, when they are long):
#
#     ```@eval
#     Main.input_file("files/bilayer/bilayer.inp")
#     ```
#
input_file(file) = Markdown.MD(Markdown.Code("", read(file, String)))
makedocs(
    modules = [Packmol],
    sitename = "Packmol.jl",
    format = Documenter.HTML(
        top_menu = true,
        # volume_tools.js: the online calculators of the Recipes tab
        assets = ["assets/custom.css", "assets/volume_tools.js"],
        # The interactive 3D views of the examples are embedded in the pages
        example_size_threshold = nothing,
        size_threshold_ignore = [
            joinpath(dir, file) for dir in ("constraints", "periodic_boundary_conditions")
            for file in readdir(joinpath(@__DIR__, "src", dir)) if endswith(file, ".md")
        ] ∪ ["tutorials/zeolite.md", "tutorials/polymer_dmf.md", "tutorials/nanotube.md", "tutorials/liquid_crystal.md", "tutorials/bilayer.md"],
    ),
    plugins = [CodeBlocks()],
    pages = [
        "Home" => [
            "Introduction" => "index.md",
            "Installation" => "installation.md",
            "Legacy (Fortran) Packmol" => "legacy.md",
        ],
        "Tutorials" => [
            "Overview" => "tutorials/index.md",
            "Using your AI agent" => "tutorials/ai_agent.md",
            "Water and urea" => "tutorials/mixture.md",
            "Water–chloroform interface" => "tutorials/interface.md",
            "Lipid bilayer" => "tutorials/bilayer.md",
            "Nanoliposomes" => "tutorials/nanoliposome.md",
            "Solvated protein" => "tutorials/solvprotein.md",
            "Protein in glycerol and water" => "tutorials/protein_glycerol.md",
            "Water drop on graphite" => "tutorials/water_drop.md",
            "Polymer in DMF" => "tutorials/polymer_dmf.md",
            "Carbon nanotube filled with water" => "tutorials/nanotube.md",
            "Gases in a zeolite" => "tutorials/zeolite.md",
            "Liquid crystals" => "tutorials/liquid_crystal.md",
            "Zika virion (a large system)" => "tutorials/zika.md",
        ],
        "Input files" => [
            "Overview" => "input_files/index.md",
            "Global keywords" => "input_files/global_keywords.md",
            "Structure blocks" => "input_files/structure_blocks.md",
            "Residue numbering" => "input_files/residue_numbering.md",
            "Restarting a run" => "input_files/restarting_a_run.md",
            "Volume estimates" => "input_files/volume_estimates.md",
        ],
        "Constraints" => [
            "Overview" => "constraints/index.md",
            "Fixed molecules" => "constraints/fixed.md",
            "Cube" => "constraints/cube.md",
            "Box" => "constraints/box.md",
            "Sphere" => "constraints/sphere.md",
            "Plane" => "constraints/plane.md",
            "Cylinder" => "constraints/cylinder.md",
            "Ellipsoid" => "constraints/ellipsoid.md",
            "Gaussian, wave, and exponential surfaces" => "constraints/gaussian_wave_and_exponential_surfaces.md",
            "Rotation constraint" => "constraints/rotation_constraint.md",
            "Julia constructors" => "constraints/julia_constructors.md",
        ],
        "Periodic boundary conditions" => [
            "Overview" => "periodic_boundary_conditions/index.md",
            "Orthorhombic box" => "periodic_boundary_conditions/orthorhombic_box.md",
            "Triclinic box" => "periodic_boundary_conditions/triclinic_box.md",
            "Rhombic dodecahedral box" => "periodic_boundary_conditions/rhombic_dodecahedral_box.md",
            "Truncated octahedral box" => "periodic_boundary_conditions/truncated_octahedral_box.md",
            "Relationship with constraints" => "periodic_boundary_conditions/relationship_with_constraints.md",
        ],
        "Julia API" => [
            "Overview" => "julia_api/index.md",
            "Multiple structure types" => "julia_api/multiple_structure_types.md",
            "Combining constraints" => "julia_api/combining_constraints.md",
            "Constraints on some atoms" => "julia_api/constraints_on_some_atoms.md",
            "Global options" => "julia_api/global_options.md",
            "Rotation constraints" => "julia_api/rotation_constraints.md",
            "Restart files" => "julia_api/restart_files.md",
            "Reference" => "julia_api/reference.md",
        ],
        "Recipes" => [
            "Overview" => "recipes/index.md",
            "Volume guesser" => "recipes/volume_guesser.md",
            "Free volume" => "recipes/free_volume.md",
            "Running Packmol" => "recipes/running_packmol.md",
            "Solute-Solvent system" => "recipes/solute_solvent_system.md",
            "Solute-Solvent-Cossolvent system" => "recipes/solute_solvent_cossolvent_system.md",
            "Solute-Water-Ions system" => "recipes/solute_water_ions_system.md",
            "Membrane system" => "recipes/membrane_system.md",
            "Membrane protein system" => "recipes/membrane_protein_system.md",
            "Concentration units" => "recipes/concentration_units.md",
        ],
    ],
)
deploydocs(
    repo = "github.com/m3g/Packmol.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
    # Pull requests from branches of this repository are published to previews/PR<number>
    push_preview = true,
)
