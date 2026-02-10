# Packmol.jl Development Guide

## Project Overview

Packmol.jl is a Julia reimplementation of the Fortran [Packmol](https://github.com/m3g/packmol/tree/master/src)
package for generating initial configurations for molecular dynamics simulations.
The goal is to achieve feature parity with the Fortran version while enabling new features,
better maintainability, and leveraging the Julia ecosystem.

The package also includes:
- A **legacy runner** (`run_packmol`) that wraps the Fortran binary via `Packmol_jll`.
- **PackmolInputCreator**: a Julia-specific module for generating Packmol input files
  for common solution setups (solute+solvent, cossolvent, ions).

## Project Structure

```
src/
  Packmol.jl                          # Main module
  packmol_runner.jl                   # Legacy Fortran binary runner
  CLI.jl                              # Command-line interface (Julia 1.12+)
  random.jl                           # Random number generation
  mono_atomic.jl                      # Monoatomic packing (simplified case)
  interatomic_distance_fg.jl          # Distance function & gradient computation
  data_structures/
    atoms_and_molecules.jl            # AtomData, MoleculePosition, FixedMoleculeData
    StructureType.jl                  # Structure definitions & input parsing
    PackmolSystem.jl                  # Main system container
  constraints/
    constraints_base.jl               # Abstract Constraint type & interface
    boxes.jl                          # Cube, Box (Inside/Outside)
    spheres.jl                        # Sphere (Inside/Outside)
  rigid_body/
    rigid_body.jl                     # Euler angles, rotations, translations
    chain_rule.jl                     # Gradient chain rule for rigid bodies
  PackmolInputCreator/                # Input generation for common solution types
    PackmolInputCreator.jl
    SolutionBoxUS.jl                  # Solute + Solvent
    SolutionBoxUSC.jl                 # Solute + Solvent + Cossolvent
    SolutionBoxUWI.jl                 # Solute + Water + Ions (incomplete)
    DensityTable.jl
    concentration_units.jl
```

## Key Dependencies

- `CellListMap.jl`: Efficient pairwise distance computation with cell lists
- `SPGBox.jl`: Spectral projected gradient optimizer (replaces GENCAN from Fortran)
- `PDBTools.jl`: PDB and mmCIF (which will be a new feature) file reading/writing
- `StaticArrays.jl`: Fixed-size vectors/matrices for performance
- `Unitful.jl`: Physical unit handling (in PackmolInputCreator)

## Running Tests

```bash
julia --project -e 'using TestItemRunner; @run_package_tests'
```

Individual test items are defined inline with `@testitem` in the source files.

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
- [~] Gradient chain rule (Cartesian -> rigid-body DOFs) - implemented, recently fixed, **needs tests**
- [ ] 2D rigid body chain rule tests

#### Input Parser
- [x] Global keywords: `tolerance`, `output`, `filetype`, `seed`, `maxit`, `writeout`, `writebad`,
      `randominitialpoint`, `avoid_overlap`, `add_amber_ter`, `amber_ter_preserve`, `add_box_sides`,
      `connect`, `optimization_print_level`, `chkgrad`, `tolerance_precision`, `constraint_precision`
- [x] Structure blocks: `structure/end structure`, `number`, `fixed`, `center`, `radius`
- [x] Constraint parsing: `inside/outside box`, `inside/outside cube`, `inside/outside sphere`
- [x] Per-atom blocks (`atoms ... end atoms`) with custom constraints
- [ ] Missing keyword parsing: `discale`, `movefrac`, `movebadrandom`, `sidemax` (won't be supported),
      `precision`, `fbins` (won´t be supported), `fscale`, `use_short_tol`, `short_tol_dist`, `short_tol_scale`, `short_radius`, `short_radius_scale`, `resnumbers`, `changechains`, `chain`, `restart_from`, `restart_to`, `constrain_rotation`, `nloop`

---

### Phase 2: Constraint System (Geometric Regions)

Each constraint type needs: data structure, penalty function, gradient, parsing, and tests.

#### Implemented
- [x] `Box` (InsideBox / OutsideBox) - axis-aligned rectangular box
- [x] `Cube` (InsideCube / OutsideCube) - axis-aligned cube
- [x] `Sphere` (InsideSphere / OutsideSphere) - sphere

#### Missing Constraints
- [ ] `Cylinder` (Inside / Outside) - center, axis direction, radius, length
- [ ] `Ellipsoid` (Inside / Outside) - center, semi-axes, orientation
- [ ] `Plane` (Over / Below) - point + normal vector (`Over`/`Below` placement types are defined but unused)
- [ ] Combined constraints per structure (already supported by architecture, needs testing)


#### Periodic boundary conditions
- [ ] Orthorhombic boxes
- [ ] Triclinic systems
- [ ] New feature: Recipes for octahedric/icosaedric boxes.

---

### Phase 3: Objective Function & Optimization

#### Distance Computation
- [x] Pairwise distance function & gradient (`cartesian_fg!`)
- [x] CellListMap integration for efficient neighbor search
- [~] Full `fg!` function combining distance + constraint penalties - **needs thorough testing**
- [ ] Constraint penalty integration into objective function (per-atom constraint evaluation)

#### Optimization Loop
- [~] Basic SPGBox-based optimization (`packmol()` function exists in `interatomic_distance_fg.jl`)
- [ ] Multi-start strategy: initial placement + iterative repositioning of bad molecules
- [ ] `movebadrandom` heuristic: randomly reposition molecules with worst overlaps
- [ ] `movefrac` control: fraction of bad molecules moved per iteration
- [ ] Convergence criteria based on `tolerance_precision` and `constraint_precision`
- [ ] Iteration output: write intermediate configurations every `writeout` steps
- [ ] `writebad`: write output even when packing fails

**Note on optimizer choice**: The Fortran Packmol uses GENCAN (large-scale bound-constrained
optimizer with truncated Newton steps). The Julia version uses SPGBox (spectral projected gradient).
We have full control over the SPGBox.jl package, if some tuning is required. But given the approximate nature of the optimization required, it is probably good enough. Otherwise we can test any of the many other solvers available in Julia, but this is for the future.

---

### Phase 4: Output Generation

#### PDB Output
- [ ] Write packed coordinates to PDB file
- [ ] Residue numbering schemes (`resnumbers` 0/1/2/3)
- [ ] Chain identifier control (`changechains`, `chain`)
- [ ] AMBER TER records (`add_amber_ter`, `amber_ter_preserve`)
- [ ] CONECT record preservation (`connect`) - requires PDBTools.jl to support connectivity.
- [ ] CRYST1 record with box dimensions (`add_box_sides`)

#### Other Output Formats
- [ ] Tinker XYZ format (won't be supported)
- [ ] XYZ format
- [ ] Moldy format (won't be supported)

---

### Phase 5: Advanced Features

- [ ] Restart capability (`restart_from`, `restart_to`)
- [ ] Short-range tolerance (`use_short_tol`, `short_tol_dist`, etc.)
- [ ] Rotation constraints (`constrain_rotation`)
- [ ] Gradient checking mode (`chkgrad` - compare analytical vs numerical gradients) - This is not necessary in the code itself, but for testing, and if the user implements new constraints, we will provide the directives on how to test their gradients. 
- [ ] 2D packing support (data structures support D=2, needs testing and validation)

---

### Phase 6: End-to-End Integration & Testing

#### First Running Example
- [ ] Pack water molecules in a box (simplest realistic test case)
- [ ] Compare output with Fortran Packmol results: not needed. The Fortran code has already a test suite based on CellListMap and Julia. We should use something similar. 
- [ ] Verify minimum distances satisfy tolerance - and the constraints. 

#### Integration Tests
- [ ] Full pipeline: parse input -> setup system -> optimize -> write output
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

### Phase 7: PackmolInputCreator (Julia-specific)

- [x] `SolutionBoxUS` - solute + solvent setup
- [x] `SolutionBoxUSC` - solute + solvent + cossolvent setup
- [~] `SolutionBoxUWI` - solute + water + ions setup (incomplete, `write_packmol_input` has early return)
- [x] Concentration unit conversions (Molarity, Molality, MoleFraction, etc.)
- [x] DensityTable interpolation

---

### Phase 8: Polish & Release

- [ ] Resolve Project.toml merge conflict (from v0.1.13 merge)
- [ ] Documentation (Documenter.jl, hosted on GitHub Pages)
- [ ] CI/CD pipeline for tests
- [ ] Package registration
- [ ] Migration guide from Fortran Packmol

---

## Current Priority: Getting a First Running Example

The most impactful next step is achieving a complete end-to-end run. This requires:

1. **Verify chain rule gradients** (tests are empty in `chain_rule.jl`)
2. **Complete the `fg!` function** to properly combine distance + constraint penalties
3. **Implement PDB output writer** for packed coordinates
4. **Wire up the full pipeline**: `read_packmol_input` -> initialize positions -> optimize -> write output
5. **Test with a simple water box** (test data already exists in `test/data/water.pdb`)

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

## Conventions
- Dimensions are parameterized: `D=2` or `D=3` (though 3D is the primary target)
- Coordinates in Angstroms, consistent with PDB format
- Tests use `@testitem` from TestItems.jl for inline test definitions
