# Migrating from legacy Packmol

Packmol.jl reads the input files of the legacy (Fortran) Packmol, and the `packmol`
command runs them in the same way:

```bash
packmol < input.inp
packmol input.inp
packmol -i input.inp -o output.pdb
```

Most input files run unchanged. This page lists the differences: what changed, what is no
longer supported, and what is new. The legacy program is still available (see
[Running the legacy Packmol](@ref)).

## What stays the same

- The syntax of the input files: the global keywords, the `structure ... end structure`
  blocks, the constraints (`inside box`, `outside sphere`, `above plane`, ...), the
  `atoms ... end atoms` blocks, `fixed`, `center`, `radius`, `resnumbers`, `chain`,
  `segid`, `add_amber_ter`, `amber_ter_preserve`, `connect`, `add_box_sides`, ...
- The conventions of the constraints: `inside cube x y z side` from the lower corner,
  `pbc a b c` as the box from the origin to `(a, b, c)`, the Euler angles of `fixed`.
- The file formats PDB and XYZ (`filetype`), the output file `<output>_FORCED` (a solution
  that satisfies only the constraints) when the packing does not converge, and the exit
  codes: `0` success, `173` not converged, `171` error in the input, `172` a file could not
  be opened, `170` other errors (and `174` for errors in the command line).

## What changed

- **Paths are relative to the input file.** The structure files, the `output` and the
  `restart_from` file of an input file are found relative to the directory of the input
  file, not to the current directory. With the input redirected (`packmol < input.inp`)
  there is no input file, and the paths are relative to the current directory, as in
  the legacy Packmol.
- **Restarts.** `restart_from` reads a coordinate file (PDB, XYZ or mmCIF, for example the
  output of a previous run), not the restart files of the legacy Packmol (centers of mass
  and Euler angles), which are not supported. `restart_to` is ignored, with a warning: the
  output file is the restart file. See [Restarting a run](@ref).
- **`constrain_rotation`** keeps the angle between an axis (`x`, `y` or `z`) and a vector
  of the molecule within a range: the longest axis of the molecule, or the vector between
  the centers of two sets of atoms. In the legacy Packmol, it bounded the Euler angles of
  the molecule, whose meaning is hard to predict. Input files that use it must be
  revised. See [Rotation constraint](@ref).
- **`centerofmass`** (of fixed molecules) uses the masses of the atoms. The legacy Packmol
  used unit masses, that is, the geometric center (the same as `center`).
- **Convergence.** The packing has converged when no two atoms of different molecules
  overlap by more than `tolerance_precision` (0.1 Å by default), with the overlap measured
  with the radii of the atoms (`tolerance/2` by default, or `radius`), and the constraints
  are satisfied within `constraint_precision`.
- **Misspelled keywords** in structure blocks are errors (the legacy Packmol ignored them,
  which could pack zero molecules of a type without notice).
- **The configurations differ.** The algorithm (the optimizer, the initial placement and
  the random number generator) is new, so the same input and seed give a different, but
  equally valid, configuration. Packmol.jl is also faster, and multithreaded (`-t N`
  limits the number of threads).

## What is no longer supported

These keywords are accepted, with a warning, and ignored, so that old input files still
run:

- `fscale`, `fbins`, `precision`, `sidemax`, `maxmove`, `nloop0`, `packall`,
  `disable_movebad`, `movebadrandom`, `iprint1`, `iprint2`, `iprint3`;
- `use_short_tol`, `short_tol_dist`, `short_tol_scale`, `short_radius`,
  `short_radius_scale`;
- `randominitialpoint`, `optimization_print_level`, `chkgrad`;
- `writebad`: the best solution is written when the packing does not converge;
- `restart_to`, `writecrd` (CHARMM CRD output), `non_standard_conect`,
  `hexadecimal_indices` (atom indices above 99999 are always written in hexadecimal);
- `nloop` inside a structure block (use the global `nloop`).

`check` still works, but `check_initial_placement` is its new name. Not supported (an
error): the formats `filetype tinker` and `filetype moldy`, and the constraints
`over xygauss`/`below xygauss`, replaced by the more general `above gaussian`/`below
gaussian` (see [Gaussian, wave, and exponential surfaces](@ref)).

## What is new

In input files and on the command line:

- mmCIF files (`.cif`), for structures, outputs and restarts;
- `-o` (the output file) and `-s` (the seed), to build replicas of a system from one
  input file (see [Replicas of a system](@ref));
- `check_initial_placement`, to see where the constraints put the molecules before packing
  (see [Checking the initial placement](@ref check_initial_placement));
- periodic cells: `pbc xmin ymin zmin xmax ymax zmax`, triclinic cells (`unitcell a b c α β γ`),
  rhombic dodecahedra and truncated octahedra (`pbc dodecahedral`, `pbc octahedral`);
  see [Periodic boundary conditions](periodic_boundary_conditions/index.md);
- new constraints: curved surfaces (Gaussian, waves and exponentials), and the new
  `constrain_rotation`;
- several conformations of a molecule, from a structure file with several models;
- `--guess-volume` and `--free-volume`, to estimate the volumes of the system and of the
  regions of its molecules (see [Volume estimates](@ref)).

From Julia only: systems defined without input files (see
[Defining systems in Julia](julia_api/index.md)), the [Recipes](recipes/index.md) (solvated
systems and membranes from target densities and concentrations), new constraints
defined by the user ([Custom constraints](constraints/custom.md)), and the packed atoms in memory (`get_atoms`), for analysis with
[PDBTools.jl](https://github.com/m3g/PDBTools.jl) and other packages.

## Running the legacy Packmol

The legacy Packmol, the Fortran program, is still available from
[m3g/packmol](http://github.com/m3g/packmol). It can also be run from Julia, with
[`run_packmol_legacy`](@ref), which works as `packmol`. It requires the `Packmol_jll`
package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
```
