# Using your AI agent

A coding agent (Claude Code, or another assistant that can write and run code on your
computer) can build a Packmol system from a description in plain language. It works best
when it writes a Julia script that uses the [Julia API](../julia_api/index.md) of Packmol.jl,
rather than an input file: in a single script it can download the structures of the
molecules (with [PDBTools.jl](https://github.com/m3g/PDBTools.jl)), select, center and
orient them, compute their masses and charges, the numbers of molecules from densities or
concentrations, and the size of the box, pack the system, check the result, and write it.
The script is a record of how the system was built, and it can be changed and run again.

This page has guidelines for the prompt, an example, and, at the end, a section written for
the agent: a summary of what Packmol.jl provides and how to use it. Point your agent to
that section in the prompt.

## Writing the prompt

Describe the system as you would to a colleague, and say what the agent should find out by
itself:

- **The molecules.** Give PDB IDs, file names, or names of molecules the agent can build
  (water, ions, common solvents); say which chains, ligands or waters of a PDB structure to
  keep. For lipids, name them: the [Membrane recipe](@ref "Membrane system") has
  conformations of POPC, POPE, POPS, POPG, DPPC, DOPC and DMPC.
- **The amounts.** Numbers of molecules, or what determines them: a density, a
  concentration (mol/L, molality, mole fraction), a salt concentration, the neutralization
  of a charged protein, an area per lipid. The agent can compute the numbers.
- **The box.** Its shape (cubic, orthorhombic, rhombic dodecahedron, truncated
  octahedron), and its size, or a margin around a solute. Periodic boundary conditions are
  the usual choice for a simulation.
- **The regions and orientations.** Where each kind of molecule goes (inside a sphere,
  above a plane, in a layer), and how molecules are oriented (lipid heads toward the water,
  molecules aligned along an axis).
- **The output.** The file and its format (PDB or XYZ), and what you will do with
  it (a molecular dynamics simulation with a given force field, for example).
- **The workflow.** Ask the agent to check the initial placement before packing, to verify
  the result (no overlaps, the molecules in their regions, the densities), and to keep the
  script.

If something is not clear in your request (the protonation of a protein, how much margin
around it), say how you want it decided, or ask the agent to ask you before deciding.

## An example prompt

```text
Write a Julia script that uses Packmol.jl to build a system for a molecular
dynamics simulation of lysozyme (PDB 1AKI, without crystal waters) in water
with 0.15 M NaCl, neutralizing the charge of the protein, in a periodic rhombic
dodecahedron with at least 12 Å between the protein and the faces of the cell.
Use PDBTools.jl to download the structure and to compute what you need. Check
the initial placement first, then pack, and verify that the water has a
density close to 1 g/mL. Write the system to lysozyme.pdb, in the shape of the
dodecahedron. Read the section "For AI agents" of this page before starting:
<paste here the address of that section>
```

## For AI agents

This section is written for an AI agent that builds systems with Packmol.jl. It summarizes
the package and how to use it in practice; the rest of this documentation has the details.

### Approach

- Write a Julia script, and run it with `julia --project=<environment> -t auto script.jl`
  (`-t` sets the number of threads; packing is multithreaded). Install with
  `import Pkg; Pkg.add(["Packmol", "PDBTools"])`.
- Build the system with the Julia API (`structure_type`, `PackmolSystem`, `packmol`) or a
  recipe, not by writing an input file, unless an input file is asked for (then
  `write_packmol_input` writes one from a recipe, or write the file directly; the syntax is
  that of Fortran Packmol, see [Input files](../input_files/index.md)).
- Use PDBTools.jl for the structures: `wget("1AKI", "protein")` downloads a structure
  (with a selection), `read_pdb`/`read_mmcif` read files, `write_pdb` writes them, `mass`
  gives masses, `center_of_mass`, `maxmin` the extent of a structure, `moveto!` moves it,
  and `select` selects atoms. `structure_type` takes the name of a structure file, so write
  each molecule to a file first.
- Units: Å for lengths and coordinates, g/mol for masses, g/mL for densities.

### Core API

```julia
using Packmol
water = structure_type("water.pdb";
    number = 1000,                                     # number of molecules
    constraints = [InsideBox([0., 0., 0.], [40., 40., 40.])],    # whole molecule
    atom_constraints = [[1] => AbovePlane([0., 0., 1.], 10.)],  # some atoms only
    constrain_rotation = Dict(:z => (0.0, 20.0, [1], [2])),     # orientation (atoms 1 to 2)
)
protein = structure_type("protein.pdb"; number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),  # position and Euler angles (radians)
    center = :geometric,                   # or :mass, or :none (coordinates as given)
)
system = PackmolSystem([protein, water];
    output = "system.pdb",          # PDB, or XYZ (by the extension, .xyz)
    tolerance = 2.0,                # minimum distance between atoms of different molecules (Å)
    unitcell = [50. 0. 0.; 0. 50. 0.; 0. 0. 50.],   # periodic cell (columns: cell vectors), centered at the origin
    check_initial_placement = true, # write the initial placement only, without packing
)
packmol(system)       # packs, writes `output`, returns `system`
system.status         # :packing_ready, :failed_packing, :failed_initial_approximation, :not_packed
atoms = get_atoms(system)   # the packed atoms (PDBTools atoms), for analysis or visualization
```

`packmol("input.inp")` packs an input file. `packmol(system; nloop, maxit, nthreads)`
changes the number of loops, the iterations per loop and the threads.

### Constraints

Of whole molecules (`constraints`) or of some atoms of each molecule (`atom_constraints`,
pairs of atom indices and constraints); the input-file keywords are in parentheses:

- `InsideBox(lo, hi)`, `OutsideBox` (`inside box xmin ymin zmin xmax ymax zmax`)
- `InsideCube(corner, side)`, `OutsideCube` (`inside cube xmin ymin zmin side`)
- `InsideSphere(center, radius)`, `OutsideSphere` (`inside sphere x y z r`)
- `InsideCylinder(center, axis, radius, length)`, `OutsideCylinder`: `center` is the
  center of one cap (`inside cylinder x y z ax ay az r l`)
- `InsideEllipsoid(center, a, b, c, scale)`, `OutsideEllipsoid` (`inside ellipsoid x y z a b c s`)
- `AbovePlane(normal, d)`, `BelowPlane`: `normal ⋅ x ≥ d` (`above plane nx ny nz d`, also `over`)
- Curved surfaces, Packmol.jl only: `AboveGaussian`, `AboveSin`, `AboveExponential`, their
  `Below...`, and the radial ones (`AboveRadialGaussian`, ...)
- `CustomConstraint(f; gradient)`, Julia only: the region `f(x) ≤ 0` of any function of
  the position (Å); the gradient by finite differences if not given, or, for example,
  `gradient = x -> ForwardDiff.gradient(f, x)`
- `constrain_rotation`: the angle between an axis (`:x`, `:y`, `:z`) and a vector of the
  molecule within `center ± halfwidth` degrees; the vector joins the centers of two sets of
  atoms, `(center, halfwidth, from, to)`, or is the longest axis of the molecule,
  `(center, halfwidth)`. This is not Fortran Packmol's `constrain_rotation`, which bounds
  Euler angles.

A structure file with several models (or a vector of files) is packed as several
conformations of the molecule, with `number` divided among them.

### Periodic cells

- `unitcell` (a 3×3 matrix, columns the cell vectors, centered at `unitcell_center`, the
  origin by default); in input files, `pbc xmin ymin zmin xmax ymax zmax` or `pbc a b c`
  (from the origin), `unitcell a b c α β γ`.
- Rhombic dodecahedron and truncated octahedron: `dodecahedral_unitcell(Float64, d)`,
  `octahedral_unitcell(Float64, d)` (`d` the distance between opposite faces; in input
  files, `pbc dodecahedral x y z d`). The output is a triclinic cell; to write it in the
  shape of the polyhedron, call `triclinic_to_dodecahedral(system)` (or
  `triclinic_to_octahedral`) before `packmol`, or set `periodic_boundary_style`.
- With periodic boundary conditions, molecules interact through the faces of the cell, and
  the constraints are evaluated at the periodic image of each molecule in the cell.

### Recipes

Systems from densities and concentrations, built directly (`packmol(recipe)`) or written as
input files (`write_packmol_input(recipe)`):

- `SolutionBoxUS`: a solute in a solvent. `SolutionBoxUSC`: a solute in a mixture of two
  solvents, at a concentration of the cossolvent (with a table of densities).
  `SolutionBoxUWI`: a solute in water with ions, at a salt concentration, neutralizing the
  charge of the solute. The shape of the cell is the `pbc` keyword of
  `write_packmol_input` and `packmol` (`:cubic`, `:orthorhombic`, `:dodecahedral` or
  `:octahedral`).
- `Membrane`: a lipid bilayer or monolayer, in water by default. Built-in lipids by name:
  `Membrane(lipids = ["POPC", "POPE"], lipid_molar_ratio = [3, 1], total_lipids = 200)`.
  Other lipids from files, with `lipid_head`, `lipid_tail` (the chain ends) and
  `area_per_lipid`.

### Volumes and numbers of molecules

- `guess_volume(input_or_system; density)`, or with pairs `"file.pdb" => number`: the
  volume of the molecules at a density, and the side of a cube.
- `free_volume(input_or_system)`: the volume of the region of each structure type, given its
  constraints and the fixed molecules, by Monte Carlo, and the density at which its molecules
  fill it. Use it to check the numbers of molecules: water fills its region at about
  1 g/mL; much more means too many molecules.
- From the command line: `packmol --guess-volume --density d input.inp` and
  `packmol --free-volume input.inp`.

### Workflow and checks

1. Compute the numbers of molecules from the volumes of their regions (subtracting solutes)
   and the densities, or use a recipe.
2. Pack with `check_initial_placement = true`: it takes a fraction of the time of the packing
   and shows whether the constraints put the molecules where intended (molecules may still
   overlap). `:failed_initial_approximation` means the constraints are incompatible or too
   tight for the molecules.
3. Pack. `:packing_ready` means that no atoms of different molecules are closer than the
   tolerance, and the constraints are satisfied. `:failed_packing` writes the best solution;
   usually there are too many molecules, or regions too tight.
4. Verify: the minimum distance between molecules, the molecules in their regions (for
   example, with PDBTools `select`), the densities, the total charge.

### Practical advice

- The tolerance of 2.0 Å is the usual one for all-atom systems.
- Rigid molecules are packed as given: prefer compact, extended conformations. Lipids pack
  fastest when straight, with the chains side by side, and with the head atoms constrained to
  a narrow band (the `Membrane` recipe does this). Bilayers of extended lipids are thicker
  than relaxed ones; an equilibration of the simulation brings them to their thickness.
- Do not constrain the same molecules with redundant regions; with periodic boundary
  conditions, a constraint for the whole cell is not needed.
- Large systems (millions of atoms) take minutes to hours and several GB of memory: check
  the initial placement first.
- The packed system has no overlaps, but it is not equilibrated: minimize its energy and
  equilibrate it at constant pressure (NPT; semi-isotropic for membranes) before a
  production simulation, with position restraints on solutes such as proteins, released
  gradually.
- `MembraneProtein` inserts a membrane protein, oriented as in the OPM database (normal
  along z, middle of the bilayer at z = 0, without its `DUM` atoms), in a bilayer of a given
  composition, with water and ions.
