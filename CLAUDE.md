# Packmol.jl Development Guide

## Project Overview

Packmol.jl is a Julia reimplementation of the Fortran [Packmol](https://github.com/m3g/packmol/tree/master/src)
package for generating initial configurations for molecular dynamics simulations.
The goal is to achieve feature parity with the Fortran version while enabling new features,
better maintainability, and leveraging the Julia ecosystem.

The package also includes:
- `packmol(input_file)` / `packmol()` (file dialog): the simple user-facing entry point, running
  the native engine. All paths in an input file (structures, `output`, `restart_from`/`restart_to`)
  are relative to its directory. (`run_packmol`, which used to run the Fortran binary, was removed.)
- A **legacy runner** (`run_packmol_legacy`) that wraps the Fortran binary via `Packmol_jll`, a weak
  dependency: the methods live in the `PackmolJLLExt` extension (`ext/`), so it requires
  `using Packmol_jll`. Packmol_jll will be dropped once the native engine replaces it.
- A **command line app** (`packmol`, via `Pkg.Apps`; `src/CLI.jl`) running the native engine.
  Not bundling a GUI toolkit is deliberate (it would dominate the size of a `juliac`-built
  binary): the file dialog (`pick_file`) shells out to OS tools instead.
- **Recipes**: higher-level, parameter-driven system setups (target densities/
  concentrations instead of molecule counts) — solute+solvent, cossolvent, water+ions,
  and membranes today, each solvation recipe with a `pbc` keyword
  (`:cubic`/`:orthorhombic`/`:dodecahedral`/`:octahedral`) selecting the periodic cell
  shape (`Membrane` is always orthorhombic); vesicles and nanotubes are planned. Each
  recipe supports both `write_packmol_input` (generate a `.inp` file) and
  `packmol(recipe; ...)` (build and pack directly, no file left behind).

## Drop-in replacement project (ongoing)

Replacing the Fortran `packmol` executable with a `juliac`-compiled Packmol.jl (feature parity,
CI builds for Linux/macOS/Windows, pip/uvx and standalone distribution) is tracked in
`dev/drop-in-replacement.md`. Read it before working on the CLI, the build or packaging, and
update its task list and session log when you stop.

## Project Structure

```
src/
  Packmol.jl                          # Main module
  packmol_runner.jl                   # run_packmol_legacy stub (methods in ext/)
  CLI.jl                              # Command-line interface, native engine
  threads.jl                          # `@_sync`: trim-safe `@sync`
  io/
    terminal.jl                       # Trim-safe stdout/stderr printing (`_print`, `_println`, `@_printf`)
    progress_bar.jl                   # Minimal progress bar (replaces ProgressMeter.jl)
    pick_file.jl                      # File dialog via OS tools (osascript/PowerShell/zenity/kdialog), run with libc `popen`
    write_output.jl                   # Output writing
  random.jl                           # Random number generation
  interatomic_distance_fg.jl          # Distance function & gradient computation
  data_structures/
    atoms_and_molecules.jl            # AtomData, MoleculePosition, FixedMoleculeData
    StructureType.jl                  # Structure definitions & input parsing
    PackmolSystem.jl                  # Main system container
  constraints/
    constraints_base.jl               # Abstract Constraint type & interface
    boxes.jl                          # Cube, Box (Inside/Outside)
    spheres.jl                        # Sphere (Inside/Outside)
    planes.jl                         # Plane (Above/Over/Below)
    cylinders.jl                      # Cylinder (Inside/Outside)
    ellipsoids.jl                     # Ellipsoid (Inside/Outside)
  rigid_body/
    rigid_body.jl                     # Euler angles, rotations, translations
    chain_rule.jl                     # Gradient chain rule for rigid bodies
  Recipes/                            # Parameter-driven system setups (recipes)
    Recipes.jl
    SolutionBoxUS.jl                  # Solute + Solvent
    SolutionBoxUSC.jl                 # Solute + Solvent + Cossolvent
    SolutionBoxUWI.jl                 # Solute + Water + Ions
    Membrane.jl                       # Lipid bilayer/monolayer
    DensityTable.jl
    concentration_units.jl
ext/
  PackmolJLLExt.jl                    # run_packmol_legacy methods (loaded with Packmol_jll)
```

## Key Dependencies

- `CellListMap.jl`: Efficient pairwise distance computation with cell lists
- `SPGBox.jl`: Spectral projected gradient optimizer (replaces GENCAN from Fortran)
- `PDBToolsBase.jl`: the core data structures of PDBTools.jl (`Atom`, ...), PDB and mmCIF
  reading/writing, and the shared function names of the M3G package ecosystem
  (`get_atoms`, ...). It is compatible with `juliac --trim`, and its types are the ones
  of PDBTools.jl (≥ 4), which re-exports them. Packmol `import`s and extends `get_atoms`
  from it rather than defining its own, so the other packages of the ecosystem
  (MolSimToolkit.jl, ComplexMixtures.jl, ...) can call `get_atoms` on a `PackmolSystem`
  too, and vice versa
- `StaticArrays.jl`: Fixed-size vectors/matrices for performance
- `Unitful.jl`: Physical unit handling (in Recipes)

## Running Tests

```bash
julia --project -e 'using TestItemRunner; @run_package_tests'
```

Individual test items are defined with `@testitem` in `test/`, mirroring the `src/` layout
(e.g. tests for `src/constraints/boxes.jl` live in `test/constraints/boxes.jl`). The exception
is `src/Recipes/test/`, which holds data files (and doctested paths) for the Recipes docs
examples, not test code.

---

## Development Roadmap

### Legend
- [x] Implemented and tested
- [~] Partially implemented or untested
- [ ] Not yet implemented
- [N] Not required/won´t be ported. 

---

### Phase 1: Core Infrastructure (Foundation)

These are the building blocks everything else depends on.

#### Data Structures
- [x] `AtomData` - per-atom data (molecule index, radii, constraints)
- [x] `MoleculePosition` - center of mass + Euler angles (optimization variables)
- [x] `FixedMoleculeData` - fixed molecule position/rotation
- [x] `StructureType` - molecule type definition (atoms, constraints, reference coords)
- [x] `PackmolSystem` - main system container with all settings

#### Rigid Body Mechanics
- [x] Euler angle rotation matrix (`eulermat`)
- [x] Molecule move/rotate (`move!`, `rotate!`)
- [x] Random molecule placement (`random_move!`)
- [x] Reference coordinate alignment via inertia tensor (`set_reference_coordinates!`)
- [x] Gradient chain rule (Cartesian -> rigid-body DOFs) - implemented and tested
- [ ] 2D rigid body chain rule tests

#### Input Parser
- [x] Global keywords: `tolerance`, `output`, `filetype`, `seed`, `maxit`, `writeout`, `writebad`,
      `randominitialpoint`, `avoid_overlap`, `add_amber_ter`, `amber_ter_preserve`, `add_box_sides`,
      `connect`, `optimization_print_level`, `chkgrad`, `tolerance_precision`, `constraint_precision`,
      `stall_tolerance`
- [x] Structure blocks: `structure/end structure`, `number`, `fixed`, `center`, `radius`
- [x] Constraint parsing: `inside/outside box`, `inside/outside cube`, `inside/outside sphere`,
      `inside/outside cylinder`, `inside/outside ellipsoid`, `above/over/below plane`
- [x] Per-atom blocks (`atoms ... end atoms`) with custom constraints
- [x] `radscale` (with legacy `discale` mapped to it)
- [x] Julia-native construction API (no input file needed): `structure_type(pdbfile; number, constraints, ...)`
      + `PackmolSystem(structure_types; output, ...)`, documented in `docs/src/experimental/julia_api.md`
- [x] `nloop`, `maxit` (stored in `PackmolSystem.nloop`/`max_iter`; `packmol(sys; nloop, maxit)` takes precedence)
- [x] Inline `#` comments stripped (as Fortran)
- [x] Deprecated Fortran keywords (`fscale`, `fbins`, `precision`, `sidemax`, `maxmove`, `nloop0`,
      `disable_movebad`, `writecrd`, `use_short_tol`/`short_*`, `iprint*`) accepted with a warning
- [x] `movefrac` (`packmol(sys; movefrac)` takes precedence), `movebadrandom` (no-op: bad molecules
      are always moved randomly), `segid` (written to PDB columns 73-76)
- [ ] CONECT options (`ignore_conect`, `non_standard_conect`, `hexadecimal_indices`): warned, ignored

---

### Phase 2: Constraint System (Geometric Regions)

Each constraint type needs: data structure, penalty function, gradient, parsing, and tests.

#### Implemented
- [x] `Box` (InsideBox / OutsideBox) - axis-aligned rectangular box
- [x] `Cube` (InsideCube / OutsideCube) - axis-aligned cube
- [x] `Sphere` (InsideSphere / OutsideSphere) - sphere
- [x] `Plane` (AbovePlane / BelowPlane) - point + normal vector; input keywords `above`/`over` (synonyms) and `below`, matching Fortran Packmol
- [x] `Cylinder` (InsideCylinder / OutsideCylinder) - finite capped cylinder: one end center, axis direction, radius, length
- [x] `Ellipsoid` (InsideEllipsoid / OutsideEllipsoid) - axis-aligned, center + 3 semi-axes + uniform scale factor
- [x] `Gaussian`/`Wave`/`Exponential` (Above.../Below..., Packmol.jl-only, no Fortran equivalent) - curved
      half-spaces `up·x = h(along·x)` generalizing `Plane` (degenerate to it at `amplitude=0`): a Gaussian
      bump, sine/cosine wave (`cos` = `sin` with `phase+π/2`), or exponential ramp extruded along `along`;
      see `src/constraints/gaussians.jl`, `waves.jl`, `exponentials.jl`
- [x] `RadialGaussian`/`RadialWave`/`RadialExponential` (Above.../Below..., Packmol.jl-only) - radially
      symmetric counterparts (`up·x = h(r)`, `r` = radial distance from the `up`-axis line through `center`,
      as in `Cylinder`): a smooth dome, concentric ripples, or a cone-like spike (the latter two have a
      genuine kink at `r=0`, handled by zeroing the gradient's radial component there)

#### Missing Constraints
- [ ] Combined constraints per structure (already supported by architecture, needs testing)


#### Periodic boundary conditions
- [x] Orthorhombic boxes (`pbc` keyword: 3 values = side lengths centered at origin; 6 values = xmin/xmax with center at midpoint)
- [x] Triclinic systems (`unitcell` keyword with CRYST1-style a b c α β γ, centered at origin)
- [x] Reference center for PBC: constraints evaluated on atoms wrapped to cell centered at `unitcell_center`
- [x] Rhombic dodecahedron box (`pbc dodecahedral cx cy cz d` keyword, `dodecahedral_unitcell`
      Julia API helper, `triclinic_to_dodecahedral`/`dodecahedral_to_triclinic` coordinate
      conversion helpers — see `src/data_structures/dodecahedron.jl`); recipes support it via
      their `pbc = :dodecahedral` keyword.
- [x] Truncated octahedron box (`pbc octahedral cx cy cz d` keyword, `octahedral_unitcell`
      Julia API helper, `triclinic_to_octahedral`/`octahedral_to_triclinic` coordinate
      conversion helpers — see `src/data_structures/octahedron.jl`, mirroring
      `dodecahedron.jl`; both share the same "nearest periodic image" search in
      `src/data_structures/periodic_cells.jl`); recipes support it via their
      `pbc = :octahedral` keyword.
- [ ] New feature: Recipes for icosaedric boxes.

---

### Phase 3: Objective Function & Optimization

#### Distance Computation
- [x] Pairwise distance function & gradient (`cartesian_fg!`)
- [x] CellListMap integration for efficient neighbor search
- [x] Full `fg!` function combining distance + constraint penalties - tested with numerical gradient check
- [x] Constraint penalty integration into objective function (per-atom constraint evaluation)

#### Optimization Loop
- [x] SPGBox-based optimization (`packmol()` function)
- [x] Multi-start strategy: initial placement + iterative repositioning of bad molecules (`movebad!`)
- [x] `movefrac` control: fraction of bad molecules moved per iteration
- [x] `<output>_FORCED`: constraint-only adjustment of the best solution, written when packing fails
- [ ] Convergence criteria based on `tolerance_precision` and `constraint_precision`
- [ ] Iteration output: write intermediate configurations every `writeout` steps
- [ ] `writebad`: write output even when packing fails

**Note on optimizer choice**: The Fortran Packmol uses GENCAN (large-scale bound-constrained
optimizer with truncated Newton steps). The Julia version uses SPGBox (spectral projected gradient).
We have full control over the SPGBox.jl package, if some tuning is required. But given the approximate nature of the optimization required, it is probably good enough. Otherwise we can test any of the many other solvers available in Julia, but this is for the future.

---

### Phase 4: Output Generation

#### PDB Output
- [x] Write packed coordinates to PDB file
- [x] `get_atoms(::PackmolSystem)` - packed atomic coordinates as `Vector{Atom}`, generated
      fresh from `molecule_positions` on every call (never cached on the system); `write_output`
      is now a thin wrapper around it. `get_atoms` itself is `import`ed (not defined) from
      `PDBToolsBase.jl` (a shared dependency across the M3G package ecosystem — see its
      own `function get_atoms end` stub) and extended here with a `PackmolSystem` method,
      rather than declared as an unrelated function of the same name; also `export`ed from
      Packmol for direct use after `using Packmol`. Under PBC, wraps each molecule rigidly by
      its center of mass (never atom-by-atom, so a molecule is never torn across a periodic
      boundary), in the shape selected by `PackmolSystem.periodic_boundary_style` (`:triclinic`,
      the default, `:dodecahedral`, or `:octahedral` — see
      `triclinic_to_dodecahedral(::PackmolSystem)`/`dodecahedral_to_triclinic(::PackmolSystem)`
      and their `octahedral` counterparts).
- [x] Residue numbering schemes (`resnumbers` 0/1/2/3, per structure type, with Fortran's
      single-/multi-residue-template auto-detected default when unset)
- [x] Chain identifier control (`changechains`, `chain`, per structure type, with Fortran's
      auto-cycling A-Z/1-9/0 default when unset)
- [x] AMBER TER records (`add_amber_ter`, `amber_ter_preserve`)
- [ ] CONECT record preservation (`connect`) - requires PDBTools.jl to support connectivity first.
- [x] CRYST1 record with box dimensions (`add_box_sides`, or automatically when PBC is set)

#### Other Output Formats
- [ ] Tinker XYZ format (won't be supported)
- [ ] XYZ format
- [ ] Moldy format (won't be supported)

---

### Phase 5: Advanced Features

- [x] Restart capability (`restart_from`, `restart_to`, whole-system and per structure type). Raw
      Packmol-format restart files only for now (a PDB-based restart, rigid-body-aligning each
      molecule's template to observed atom positions, was discussed and judged feasible but not
      yet built).
- [N] Short-range tolerance (`use_short_tol`, `short_tol_dist`, etc.): deprecated, warned and ignored
- [x] Rotation constraints (`constrain_rotation`) - implemented as hard bounds on the rotation-angle
      optimization variables (matching Fortran's own `pgencan` bound-setting), not a soft penalty.
- [ ] Gradient checking mode (`chkgrad` - compare analytical vs numerical gradients) - This is not necessary in the code itself, but for testing, and if the user implements new constraints, we will provide the directives on how to test their gradients. 
- [ ] 2D packing support (data structures support D=2, needs testing and validation)

---

### Phase 6: End-to-End Integration & Testing

#### First Running Example
- [x] Pack water molecules in a box (simplest realistic test case)
- [ ] Compare output with Fortran Packmol results: not needed. The Fortran code has already a test suite based on CellListMap and Julia. We should use something similar. 
- [x] Verify minimum distances satisfy tolerance - and the constraints.

#### Integration Tests
- [x] Full pipeline: parse input -> setup system -> optimize -> write output
- [x] PBC packing: periodic boundary conditions with atom wrapping
- [ ] Multi-structure systems (e.g., protein + water + ions)
- [ ] Fixed molecule + mobile molecules
- [ ] Multiple constraint types in one system

#### Regression Tests Against Fortran Packmol
- [ ] Water box
- [ ] Protein solvation
- [ ] Lipid bilayer setup
- [ ] Mixed solvent systems

#### Performance Benchmarks
- [ ] Compare packing time with Fortran Packmol for standard benchmarks
- [ ] Profile and optimize hot loops
- [ ] Memory usage comparison

---

### Phase 7: Recipes (Julia-specific)

- [x] `SolutionBoxUS` - solute + solvent setup
- [x] `SolutionBoxUSC` - solute + solvent + cossolvent setup
- [x] `SolutionBoxUWI` - solute + water + ions setup (target ionic concentration, automatic solute charge neutralization)
- [x] `Membrane` - lipid bilayer/monolayer setup, one or more lipids (`lipid_molar_ratio`) solvated by
      one or more solvents (`solvent_molar_ratio`), always orthorhombic PBC with the membrane normal
      along z. Each lipid's single head/tail atom indices give its intrinsic head-to-tail length (the
      longest one sets the membrane's nominal thickness); each copy is packed with two whole-molecule
      `plane` constraints bounding it to its leaflet's slab plus two per-atom `plane` constraints (via
      `atoms <idx> ... end atoms`) pinning head/tail near their target planes within `flexibility` (a
      fraction of that lipid's own length) for tilt/wobble. Box lateral size and lipid count are tied
      together via `area_per_lipid` plus exactly one of `total_area`/`total_lipids` (mirroring
      `box_sides`/`margin` elsewhere); solvent slab thickness/density/mixture are independent knobs.
      Every parameter is given up front in the constructor (no separate box-sizing args at
      `write_packmol_input`/`packmol` call time, unlike the solvation recipes below).
- [x] `packmol(recipe; ...)` - build and pack a recipe directly, without an explicit input file
- [x] Concentration unit conversions (Molarity, Molality, MoleFraction, etc.)
- [x] DensityTable interpolation
- [x] Documented on its own `experimental/recipes.md` docs page
- [x] `pbc` keyword (`:cubic`/`:orthorhombic`/`:dodecahedral`/`:octahedral`) selecting the
      periodic cell shape, shared by the three solvation recipes; cell volume computed generally as
      `det(unitcell)` rather than assuming an orthorhombic box
- [ ] Future recipes: vesicles, nanotubes

---

### Phase 8: Polish & Release

- [x] Resolve Project.toml merge conflict
- [ ] Documentation (Documenter.jl, hosted on GitHub Pages)
- [ ] CI/CD pipeline for tests
- [ ] Package registration
- [ ] Migration guide from Fortran Packmol

---

## Architecture Notes

### Optimization Variables
Molecules are parameterized by `MoleculePosition{D,T}` containing center-of-mass (`cm`)
and Euler angles (`angles`). The vector of `MoleculePosition` is reinterpreted as a flat
`Float64` array for the optimizer via `reinterpret`.

### Constraint Interface
All constraints implement:
- `constraint_penalty(c::Constraint, x::SVector)` -> penalty value
- `constraint_gradient(c::Constraint, x::SVector)` -> gradient vector

Penalties are quadratic: `weight * (violation)^2` with default weight 5.0.

### Cell Lists
Distance computation uses `CellListMap.jl` which provides efficient O(N) pairwise
computation with automatic cell list construction. This replaces the custom cell list
implementation in the Fortran version.

## `juliac --trim` compatibility

The CLI (`src/CLI.jl`) compiles with `juliac --trim=safe` (Julia 1.13+: with 1.12 the
bodies of spawned tasks are not compiled into the executable), which requires every call
reachable from it to be resolved at compile time. In the engine code:
- print with `_print`/`_println`/`@_printf` (`src/io/terminal.jl`), never `println`/`@printf`
  to `stdout`/`stderr` (non-constant globals of abstract type `IO`, and not even initialized
  in a trimmed executable: the CLI sets `_use_core_streams[] = true`); progress with `ProgressBar`;
- use `@_sync` (`src/threads.jl`) instead of `@sync`, and no `deepcopy`;
- read/write PDB files with PDBToolsBase's `read_pdb`/`write_pdb` (trim-compatible);
- call functions on constraints taken from a `Vector{AnyConstraint{T}}` through
  `_constraint_dispatch(f, c, args...)` (the union is too large for automatic union splitting);
- don't capture in closures (incl. `do` blocks, generators, `@spawn`) variables that are
  reassigned (they get boxed); variables captured by `@spawn` must have concrete types
  (not `Union{Nothing,...}`: use a function barrier or a flag); don't call with
  `Val(runtime_value)`, keyword calls with `Union`-typed arguments, or a `Matrix` unit cell
  where the dimension must be static (use `SMatrix{D,D,T}`);
- no `readlines` (use a `readline` loop); pass a caught exception (typed `Any`) only to
  `@nospecialize` functions, and narrow `AbstractString` fields (`e.msg`) to `String` before printing.
To check (`juliac` is a `Pkg.Apps` app; `JULIA_APPS_JULIA_CMD` selects the Julia it runs with):
`JULIA_APPS_JULIA_CMD=$(which julia) juliac --output-exe packmol --project . --experimental --trim=safe --bundle build app.jl`,
with `app.jl` defining `(@main)(ARGS) = Packmol.CLI.main(ARGS)`. Requires the trim-compatible
versions of CellListMap (`show_progress = nothing` default, typed task waits) and SPGBox
(`iprint = nothing` default).

## Conventions
- Dimensions are parameterized: `D=2` or `D=3` (though 3D is the primary target)
- Coordinates in Angstroms, consistent with PDB format
- Tests use `@testitem` from TestItems.jl, placed in `test/` mirroring `src/`'s layout (not inline in `src/`)
