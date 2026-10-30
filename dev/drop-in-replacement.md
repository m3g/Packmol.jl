# Packmol.jl as a drop-in replacement for Fortran Packmol

Tracking document for the multi-session project of replacing the Fortran `packmol`
executable (github.com/m3g/packmol, v21.2.2) with a `juliac`-compiled Packmol.jl.

**For agents starting a session:** read this file first. Choose an open task (`[ ]`),
mark it `[~]` with your session date while you work on it, and add a line to the
*Session log* at the bottom when you stop: what you did, what's left, and anything
surprising. Keep task IDs stable (don't renumber them). Put decisions in *Decisions*,
not only in the log. A Fortran checkout used for comparison lives at
`~/programs/packmol` on the maintainer's machine (`src/getinp.f90` is the keyword parser,
`testing/` is the test suite).

Legend: `[x]` done · `[~]` in progress · `[ ]` open · `[N]` won't do · `[?]` needs a decision

---

## Status (2026-10-30)

- **Feature parity (A): essentially done.** The CLI contract (stdin, `-i`/`-o`, exit
  codes, `_FORCED`), XYZ, CONECT, TER, `segid`, `nloop`/`maxit`/`movefrac`, deprecation
  warnings, fixed molecules and the Fortran `pbc` convention all match Fortran, and the
  whole Fortran test suite packs (`test/fortran_test_suite.jl`). Packmol.jl is faster
  than Fortran on every test input (A5.10). Deliberate differences: D1, D11–D16.
  `xygauss`, `writeout` as a time interval and `writebad` are not needed (D17).
  Remaining: keywords without effect (A5.8). The large test (A4.3) converges in 33 min.
- **Build (B):** `juliac/build.jl` builds the executable and `juliac/test_executable.jl`
  tests it (the whole test suite, the CLI and the exit codes: A4.1). The workflow
  `BuildExecutable.yml` runs both on Linux (x86_64, arm64), macOS (arm64, x86_64) and
  Windows (B3, B4: all pass). Threads: B7 done (D18). Next: B5 (manylinux), C1/C2.
- **Distribution (C):** not started (needs B3).
- **Docs (D):** the page of differences and Julia-only features is not written.
- Merged: PR #42 (CLI, formats, keywords), PR #43 (fixed molecules, restarts from
  coordinates, principal axes, radius-aware convergence, PBC trap), PR #44 (the new
  `constrain_rotation`, D15, A5.11). Branch `build-ci`: B2, B3.

---

## Decisions

- **D1 (2026-10-29): paths in an input file are relative to the input file's directory**
  (Packmol.jl behavior; Fortran uses the current directory). This is a documented
  difference from Fortran, not a bug. With input on stdin (`packmol < input.inp`) there is
  no input file, so paths resolve against the current directory, which matches Fortran in
  that case. The Fortran test inputs (`structure ./structure_files/...`, run from `testing/`)
  have to be adapted, or run through stdin.
- **D2 (2026-10-29): in the long run, the PyPI name is `packmol`**, the one the Fortran
  wheels already use (`~/programs/packmol/pyproject.toml`, `python/`,
  `.github/workflows/build-wheels.yaml`). **Swap it in only once everything is verified to
  work.** Until then, test the wheels on TestPyPI or under a temporary name.

- **D3 (2026-10-29, confirmed): keep `packmol < input.inp`.** It is the canonical
  invocation in the Packmol user guide and in tutorials, and it is what scripts and
  third-party wrappers use. Dropping it would make this "almost drop-in". Supporting it
  costs little: read stdin only when no argument is given and stdin is not a TTY, and keep
  the file dialog for the interactive case. The docs can promote `packmol input.inp`/`-i`
  as the preferred form.

- **D4 (2026-10-29): keywords removed from Packmol warn, they never stop execution.**
  The short-tolerance options (`use_short_tol`, `short_tol_*`, `short_radius*`) are no
  longer needed, so they are deprecated.
- **D5 (2026-10-29): XYZ I/O goes into PDBToolsBase** (it's simple and belongs with the
  readers and writers there). Packmol then uses it for `filetype xyz`.
- **D6 (2026-10-29): CONECT goes through PDBToolsBase first.** PDBToolsBase reads the
  connectivity into a `Vector{Vector{Int32}}` in the new `meta_data` field of `Structure`,
  and Packmol adapts it (offsetting the indices per molecule copy) and forwards it to its
  output.
- **D7 (2026-10-29): `nloop`/`maxit` are read into the system** (`sys.nloop`,
  `sys.max_iter`). The `nloop`/`maxit` keyword arguments of `packmol(sys; ...)` take
  precedence.

- **D8 (2026-10-29): exit codes are Fortran's** (0, 170–174). They are specific to
  Packmol. The closest thing to a convention is BSD `sysexits.h` (64 usage, 65 data error,
  66 no input, 70 software error), which is no real standard either, so we keep Fortran's
  codes for drop-in compatibility.
- **D9 (2026-10-29): `segid` is written to PDB columns 73–76** for every atom of the type.
  Fortran used it only in CRD output (which won't be supported). Fortran also copied the
  `segid` of the previous structure block into a block that has none (with a warning).
  That was a bug in Fortran (maintainer, 2026-10-29). Here a block without `segid` keeps the
  segment names of its own PDB file, or blank if it has none.
- **D10 (2026-10-29): `movebadrandom` is deprecated** (warned and ignored). Julia's `movebad!` always
  places bad molecules at random positions (within the per-type center-of-mass bounds,
  with a short constraint fit), which is Fortran's `movebadrandom` behavior. Fortran's
  default instead places them next to a randomly chosen good molecule. That mode is not
  implemented; add it only if benchmarks show it helps.

- **D11 (2026-10-29): `pbc a b c` is the box from the origin to (a, b, c)**, as in
  Fortran (`pbc 0 0 0 a b c`). Packmol.jl used to center it at the origin, which shifted
  every position-dependent constraint of Fortran inputs by half a box. That's why
  `bilayer_pbc.inp` and `water_box_pbc2.inp` didn't converge. `unitcell a b c α β γ`
  (Packmol.jl only, used by the recipes) stays centered at the origin.
- **D12 (2026-10-29): the minimum-distance criterion is Packmol.jl's own**
  (`tolerance - dmin < tolerance_precision`, default 0.1 Å). It is deliberately looser
  than the 0.99 × tolerance that Fortran's `runtests.jl` checks.
  *Refined 2026-10-29 (to confirm):* the overlap of each pair is measured with the
  radii of its atoms, `max(r_i + r_j - d) < tolerance_precision`. That's the same
  thing with the default radii (tolerance/2), but `radius` larger than tolerance/2 is
  now respected (it was ignored: overlaps of 0.25 Å were accepted), and a smaller
  `radius` allows closer pairs, as in Fortran.
- **D13 (2026-10-29): no Fortran restart files.** Restarts use coordinate files only
  (PDB, XYZ or mmCIF, by extension): `restart_to` writes the coordinates of the system
  (or of one structure type), and `restart_from` recovers each molecule's position from
  coordinates by rigid-body alignment to its template. The raw format of center of mass
  and Euler angles was removed, so Fortran restart files are not read (a clear input
  error), and Fortran never reads Packmol.jl's.
- **D14 (2026-10-29): Packmol.jl's own Euler angles (`eulermat`) are the only convention.**
  Fortran's z–x–z `eulerrmat` of free molecules is not used anywhere.
- **D15 (2026-10-29, revised by the maintainer): `constrain_rotation` constrains the
  angle between an axis and a molecular vector**, regardless of what Fortran does
  (Fortran compatibility is not required). `constrain_rotation z <center> <halfwidth>`
  keeps the angle between the z axis and a vector of the molecule within
  `center ± halfwidth`. The vector joins the centers of two sets of atoms, given in a
  block (syntax approved by the maintainer):

      constrain_rotation z 0. 15.
        atoms 1 2 3
        to 25 26
      end constrain_rotation

  Without the block, the vector is the longest principal axis of the molecule. It is a
  penalty, like the other constraints (A5.11).
- **D16 (2026-10-29): reference orientations.** Free molecules are oriented internally
  along their principal axes of inertia (equal masses: the longest axis along z, the
  largest moment along x; `set_reference_coordinates!`), which defines their zero
  angles. Fixed molecules are rotated from
  the given coordinates: about the origin of the coordinates, or about the geometric
  center (`center`) or the center of mass (`centerofmass`, with the atomic masses,
  where Fortran uses unit masses) when asked.

- **D17 (2026-10-30): not ported.** `over/below xygauss` (the Gaussian surface will be an
  example of a custom constraint function in the Julia API), `writeout` as a time interval
  and `writebad` (A3.4, A5.6, A5.7).
- **D18 (2026-10-30): threads.** The executable starts all the threads of the CPUs
  available, and the number used is set by `-t N`, the input keyword `nthreads` or
  `JULIA_NUM_THREADS` (in this order of precedence). BLAS uses one thread (B7).
- **D19 (2026-10-30): the positional Julia constructors of boxes and cubes take the
  parameters of the input file** (maintainer's choice): `InsideBox(lower, upper)` and
  `InsideCube(lower_corner, side)`. The keyword forms (`center`, `sides`/`side`) keep the
  center convention. Before, the positional forms took the center, unlike the input lines
  and the docs.
- **D20 (2026-10-30): distribution (maintainer).** Everything lives in this repository,
  including the Python packaging; m3g/packmol will be archived (or renamed, with this
  repository's content merged into it, for continuity). The maintainer provides only
  GitHub release archives (C1) and PyPI/uvx wheels (C2). Linux distribution packages,
  Homebrew, conda-forge, etc. are left to third parties, who can repackage the C1
  archives. The installation page of the docs links to the binaries of the latest
  release and of the development version.
---

## Baseline (measured 2026-10-29, Linux x86_64, Julia 1.13.2)

- `juliac --trim=safe --bundle` builds the CLI successfully (with a warm depot it takes
  about 16 s and peaks at about 2 GB RAM). Command as in `CLAUDE.md`, with `app.jl` =
  `using Packmol; (@main)(ARGS) = Packmol.CLI.main(ARGS)`.
- Bundle: 117 MB uncompressed (`bin/` 23 MB, `lib/julia/` 94 MB). Compressed it is
  36 MB, both as a `.tar.gz` and as a wheel-like zip (each library file stored once, with
  no symlinks). Requires glibc ≥ 2.34. **Size is not a limit for PyPI** (100 MB per file).
- Slimming experiment (B6): the binary still runs all tests after removing libgfortran,
  libquadmath, libgomp, libmpfr and libgmp, and after replacing OpenBLAS (35 MB) with a
  15 KB stub. Result: 67 MB uncompressed, 21 MB compressed. OpenBLAS and
  libblastrampoline themselves can't be deleted: LinearAlgebra's `__init__`, pulled in by
  StaticArrays, loads them at startup. The stub has to define `isamax_64_`, `zdotc_64_`,
  `cdotc_64_` and `sdot_64_` (correct implementations), which libblastrampoline probes to
  detect the interface. Packmol makes no BLAS calls. Stub source:
  `dev/blas_stub.c`.
- All 21 inputs in Fortran's `testing/input_files` parse, once paths are adapted (D1).
  The examples parse too, except `examples/interface.inp`, which uses `filetype xyz` (F1).
- Runtime against Fortran (single thread, s): water_box 0.05 vs 0.13 · mixture 0.11 vs 1.34 ·
  solvprotein 1.31 vs 0.98 · bilayer 0.04 vs 0.07 · spherical 0.11 vs 0.14.
  The binary honors `JULIA_NUM_THREADS` (solvprotein: 2.47 s with 1 thread, 0.87 s with 4,
  on a second run).
- Requirement: `Pkg.resolve()` after the PDBToolsBase merge (stale local Manifest).

---

## Workstream A: feature parity (the drop-in gaps)

The input language is almost complete. The real gaps are the CLI contract and the
coordinate file formats.

### A1. CLI contract (essential: scripts and wrappers call packmol this way)
- [x] **A1.1** `packmol < input.inp`: read the input from stdin when no argument is given
  and stdin is not a TTY. Today, with no argument, the CLI opens the file dialog and the
  call fails with exit 1. Keep the dialog only for an interactive TTY. Paths resolve
  against the current directory (D1). `read_packmol_input` needs an `IO` method.
- [x] **A1.2** `-o output.pdb` flag (overrides `output`). Fortran accepts
  `packmol -i in.inp -o out.pdb` and the arguments in any order (`src/cli_parser.f90`).
- [x] **A1.3** Exit codes as in Fortran `src/exit_codes.f90`: 0 success, 170 general error,
  171 input error, 172 file-open error, 173 failed to converge, 174 command-line error.
  Today the CLI returns 1 for all of these.
- [x] **A1.4** `<output>_FORCED`: when packing fails, the best solution goes to the output,
  and a constraint-only adjustment of it (`adjust_constraints!`) goes to `<output>_FORCED`,
  as Fortran's `checkpoint.f90` does. (`writebad` is parsed but has no other effect.)
- [x] **A1.5** `--version` / `-h` flags (useful to package managers and to their tests).

### A2. File formats
- [x] **A2.1** `filetype xyz` (2026-10-29), with `read_xyz`/`write_xyz` from PDBToolsBase 1.1.
  The format of each file follows its extension (`.pdb`, `.xyz`, `.cif`/`.mmcif`), and
  `filetype` applies otherwise. The XYZ comment line is "Built with Packmol".
  `examples/interface.inp` runs. It also exposed a bug, now fixed: a `fixed` structure
  without `number` had 0 molecules (Fortran treats it as one molecule).
- [x] **A2.2** `filetype tinker` (and `moldy`) won't be supported: any `filetype` other
  than `pdb` or `xyz` is an input error.
- [N] **A2.3** `writecrd` (CHARMM CRD output): deprecated with a warning (D4; the
  maintainer judged CRD of no use).
- [x] **A2.4** CONECT (2026-10-29), with PDBToolsBase 1.1 `Structure.connectivity`.
  - Each template's connectivity is stored in `StructureType.connectivity`, repeated
    with offset indices for every copy (`_system_connectivity`), and written by
    `write_pdb` (or after the TER records).
  - The default is on, as in Fortran. It can be turned off for one type with
    `connect no` inside its structure block, or globally with `ignore_conect` (or
    `connect no` at the top level).
  - `hexadecimal_indices` is deprecated: PDBToolsBase always writes hexadecimal serials
    above 99999.
  - `non_standard_conect` is deprecated: CONECT parsing follows the PDB standard
    columns in PDBToolsBase.
  - Output is always written with the PDBToolsBase writers (`write_pdb`/`write_xyz`).
    When TER records are requested, a second pass over the written file inserts them
    after the right atom records (`_add_ter_records`), as the maintainer asked.
  - **Requires PDBToolsBase 1.1 to be released** (compat bumped to `1.1`).
- [x] **A2.5** AMBER TER (`add_amber_ter`: a TER after every molecule; `amber_ter_preserve`:
  the template's TER records are repeated in every copy), as in Fortran `output.f90`.
  Done 2026-10-29 (`_add_ter_records` in `src/io/write_output.jl`).

### A3. Input keywords
- [x] **A3.1** Strip inline comments (`inside box 0 0 0 10 10 10  # solvent`), as Fortran does.
  Today only lines starting with `#` are skipped.
- [x] **A3.2** `nloop` and `maxit` are read from the input (D7). Before this change, the
  input's `maxit` was parsed but never used either. Note that `maxit` now means SPGBox
  iterations per loop (Julia default 800), while in Fortran it meant GENCAN iterations
  (default 20), so old inputs that set `maxit` get fewer iterations than the Julia default.
  `nloop0` and per-structure `nloop` are deprecated (there is no per-type packing phase).
  Still to check: that `water_box_failed.inp` now fails as expected (A4.1).
- [x] **A3.3** Deprecated keywords warn instead of erroring (D4): `sidemax`, `maxmove`,
  `nloop0`, `disable_movebad`, `writecrd`, the short-tolerance options, `fscale`, `fbins`,
  `precision` and `iprint*`. (`segid` and the CONECT options were implemented later:
  A3.3b, A2.4.)
- [x] **A3.3b** `segid` (D9).
- [N] **A3.4** `over xygauss` / `below xygauss` (Fortran's Gaussian surface): not needed
  (maintainer, 2026-10-30, D17). The Gaussian surface will be an example of a custom
  constraint function in the Julia API.
- [x] **A3.5** `movefrac` is read into `sys.movefrac` (the `packmol()` keyword argument
  takes precedence). `movebadrandom` is deprecated (D10).

### A4. Validation
- [x] **A4.1** Port Fortran's `testing/` suite (runtests.jl validates outputs with
  CellListMap; test.sh, test_cli.sh, test_failed.sh) to run **against the compiled
  binary**, with the inputs adapted for D1. This is the drop-in acceptance test.
  - Status (2026-10-29): the test files are now the `packmol` set of MolSimTestFiles
    (`packmol-v1`). It holds Fortran's test inputs (adapted to D1), its examples and
    the reference outputs of its connectivity test.
  - Done: `benzene2` (bond orders preserved, which needed the PDBToolsBase fix in
    m3g/PDBToolsBase.jl#7) and `examples/interface.inp` run as tests.
  - The compiled binary converges on every input except three:
    - `bilayer_pbc.inp` and `water_box_pbc2.inp` didn't converge (50 s and 110 s).
      Fixed (D11): they now converge in under 0.1 s.
    - `protein_outside_pbc_error.inp`: fixed (2026-10-29). A fixed molecule with an
      atom outside the periodic cell is now an input error (exit 171), checked when
      the input is read and when `packmol()` starts. An atom is accepted if it lies
      in the parallelepiped of the cell vectors around the cell center, or in the
      nearest-image region (dodecahedral and octahedral cells); for orthorhombic
      cells both are the same box, as in Fortran.
      - This also exposed a bug, now fixed: `_collect_fixed_positions` applied the
        fixed rotation and translation to `reference_coordinates` a second time,
        although they already hold the placed atoms. The fixed atoms used by the
        overlap checks (initial placement, `adjust_constraints!`, `movebad!`) were
        therefore misplaced for any fixed molecule not at the origin with zero
        rotation. The packing objective itself was not affected.
  - `test/fortran_test_suite.jl` (2026-10-29) packs every input of Fortran's
    `test.sh` list and checks the intermolecular minimum distance against D12. It
    also covers `water_box_failed` (`_FORCED`) and `protein_outside_pbc_error`.
    `toy_with_fixed.inp` had been written for the old centered `pbc`; it's fixed
    in `packmol-v2`.
  - Done (2026-10-30): `juliac/test_executable.jl` runs the compiled executable on every
    input of the set and checks the outputs (minimum intermolecular distance), the CLI
    (`--version`, `-h`, stdin, `-i`, `-o`, positional input) and the exit codes (171,
    173 with `_FORCED`, 174). It runs in `BuildExecutable.yml` (B3) on all platforms.
    Not ported from `test_cli.sh`: identical outputs across invocation methods.
- [x] **A4.2** Benchmark the binary against Fortran on the test inputs: see the table in
  A5.10 (Packmol.jl is faster on all of them).
- [x] **A4.3** **The large test** (run only when the maintainer asks):
  `~/Documents/Zika/montevideo/packmol/all_from_scratch.inp`, a Zika virion (fixed
  capsid, two membranes, ~550k water and ion molecules), about 30 min. Its
  `restart_to ... .pack` line is rejected since D13, so run a scratch copy (symlinked
  `pdb_files/`, that line pointed to a `.pdb`). Record the times in the session log.

---

### A5. Audit round 2 (2026-10-29): Fortran and Packmol.jl run side by side

This round compared outputs of the two binaries on the same inputs, not just keyword lists.

Fixed (branch `fixed-molecules-compat`):
- [x] **A5.1 Fixed molecules were placed twice.** `_apply_fixed_center!` already puts the
  placed atoms in `reference_coordinates`, but initialization also set the molecule
  position to `fixed.position`. The position and rotation were therefore applied a second
  time, in both the packing objective and the output, for any fixed molecule not at the
  origin with zero rotation (e.g. `fixed 5 -3 2 0 0 0` was written 5 Å off).
  - Fixed molecules now have the identity position, and their reference coordinates are
    their atoms.
  - The output does not wrap fixed molecules (they must be inside the cell).
  - `_collect_fixed_positions` (overlap checks) uses the same positions.
- [x] **A5.2 `fixed` without `center`** rotated about the geometric center; Fortran
  rotates about the origin of the coordinates (x = R·r + t). Now the same as Fortran.
  The Euler-angle convention of `fixed` (Fortran's `eulerfixed`) already matched.
- [x] **A5.3 Restart files** used Packmol.jl's angles while Fortran's use its
  `eulerrmat` (z–x–z) convention. Superseded by D13: restart files are coordinate files.
- [x] **A5.4 Misspelled keywords in structure blocks** (`numbr 10`) were silently ignored
  (0 molecules). They are now errors. Deprecated and global keywords inside structure
  blocks are ignored with a warning (Fortran accepts them there).

- [x] **A5.5 `constrain_rotation` semantics differ.** Fortran bounds the angles of its
  z–x–z `eulerrmat` convention (x → θ, y → β, z → γ, `app/packmol.f90`), so for example
  `constrain_rotation x 0 10` alone restricts the tilt of the molecule's z axis to 10°.
  Packmol.jl bounds its own Euler angles (`eulermat`). Neither is the intended behavior:
  replaced by D15 (A5.11).
- [N] **A5.6 `over/below xygauss`** constraint (= A3.4, D17).
- [N] **A5.7 `writeout`/`writebad`.** In Fortran, `writeout` is a time interval in
  seconds (default 10) for writing the best solution so far, and `writebad` also writes
  worse points. Not needed (maintainer, 2026-10-30, D17): Packmol.jl parses both and
  keeps writing after every improved loop.
- [ ] **A5.8 Parsed but without effect:** `randominitialpoint` (Fortran skips the
  constraint fit of the initial point), `optimization_print_level`, `chkgrad`.
  Implement them or deprecate them.
- [x] **A5.9 `centerofmass` of fixed molecules:** Fortran sets all masses to 1, so its
  "center of mass" is the geometric center. Packmol.jl uses atomic masses (0.05 Å
  difference for the protein test). Kept (D16). To document (D-doc).
- [x] **A5.11 New `constrain_rotation` (D15)**, `src/constraints/rotations.jl`
  (branch `rotation-constraints`). The vector is a fixed linear combination of the atom
  positions, v = Σ wᵢ xᵢ (Σ wᵢ = 0): the difference of the centers of the two sets, or,
  by default, L·R ẑ (wᵢ = L zᵢ/Σ z², L the length of the molecule along z). The penalty
  `weight·(L·Δθ)²` is the squared displacement of the tip of the vector (Å², comparable
  to the other penalties and to `constraint_precision`); its gradient goes to the atoms
  through the weights and the existing chain rule (at constant |v|, as the molecule is
  rigid). It is evaluated wherever constraints are (objective, constraint-only fits,
  badness of `movebad!`, initial placement checks). The hard bounds on the Euler angles
  (`rotation_bounds`) were removed, so the optimizer has no angle bounds left; the
  initial angles are uniform. Julia API: `constrain_rotation = Dict(:z => (c, hw))` or
  `Dict(:z => (c, hw, from, to))`.

- [x] **A5.12 `inside/outside cube x y z side` read (x, y, z) as the center** (2026-10-30).
  Fortran (and the docs) take the lower corner (`comprest.f90`: `xmin = restpars(1)`,
  `xmax = xmin + side`). Fixed in `parse_constraint`; none of the test inputs use cubes.
- [x] **A5.10 Performance: `mixture_pbc` was trapped by a PBC inconsistency.** The
  centers of mass are bounded to the closed cell, but `wrap_to_center` maps the upper
  faces to the lower ones. A urea pushed against z = +20 was evaluated for its
  constraints at z = −20, violating `over plane 0` by the whole cell, with a gradient
  the bound didn't let it follow. Packings took 10–400 s and failed with ~12% of the
  seeds. Fixed with `_wrap_cm` (centers of mass in the closed cell are not wrapped).
  Wall-clock times (s), trimmed binary vs Fortran, 2026-10-29:

  | input | Fortran | Packmol.jl |
  |---|---|---|
  | mixture_pbc | 0.70 | 0.06 |
  | mixture | 1.33 | 0.10 |
  | water_box | 0.11 | 0.05 |
  | solvprotein | 12.7 | 0.35 |
  | solvprotein_pbc | 2.73 | 1.89 |
  | spherical | 36.7 | 4.4 |
  | spherical_pbc | 0.99 | 0.62 |
  | bilayer | 26.3 | 2.5 |
  | bilayer_pbc | 0.61 | 0.60 |

  Recheck the multi-thread slowdown seen before the fix (B7).

## Workstream B: standalone executable build (and CI)

- [x] **B1** Trim-safe build works locally on Linux (see Baseline).
- [x] **B2** Build script in the repo (2026-10-30): `juliac/` holds `app.jl`, `build.jl`
  (JuliaC's library API; `julia --project=juliac juliac/build.jl [outdir]`, default
  `juliac/build`), `test_executable.jl` and a `Project.toml` (JuliaC, Packmol from `..`,
  and the test dependencies). The Julia version is pinned in the workflow (1.13.2).
- [x] **B3** CI workflow `BuildExecutable.yml`: a matrix with ubuntu-latest,
  ubuntu-24.04-arm, macos-latest (arm64), macos-15-intel and windows-latest, each building
  the bundle, running the tests (A4.1) on the bundle, and uploading it as an
  artifact. It runs the whole matrix on PRs, pushes to main, tags and
  `workflow_dispatch` (any new non-trim-safe call breaks the build). If the matrix gets
  too slow for PRs, restrict PRs to ubuntu-latest.
- [~] **B4** **Windows and macOS spikes.** It is unverified whether `juliac --trim` works
  on them with 1.13. On Windows, look at the toolchain JuliaC needs, and at `pick_file`'s
  PowerShell path. On macOS, look at code signing/notarization (an unsigned binary is
  quarantined when downloaded through a browser, though pip/conda installs are fine).
  Record the findings here.
  - Findings (2026-10-30, PR #45): `juliac --trim=safe` builds and passes all 81
    executable checks on ubuntu-latest, ubuntu-24.04-arm, macos-latest (arm64),
    macos-15-intel and windows-latest, with no platform-specific code. Windows uses the
    MinGW that JuliaC downloads by itself. Compressed bundles: Linux x86_64 37 MB, arm64
    35 MB, macOS arm64 15 MB, x86_64 21 MB, Windows 49 MB.
  - Still open: code signing/notarization on macOS, and `pick_file` (PowerShell) on
    Windows, which can't be tested in CI.
- [ ] **B5** Linux portability: the build host's glibc sets the minimum (≥ 2.34 today).
  Checked on the CI-built executable (2026-10-30): the newest glibc symbol it needs is
  `GLIBC_2.34`, from the executable itself (the bundled Julia libraries need ≤ 2.17).
  So it runs on Ubuntu ≥ 22.04, Debian ≥ 12, RHEL/Rocky/Alma ≥ 9, but not on RHEL 8
  and its clones (glibc 2.28, supported until 2029 and common on HPC clusters), Ubuntu
  20.04 or Debian 11 (2.31).
  Build in a manylinux_2_28 container (or older) for wider coverage, matching the wheel tag.
- [~] **B6** Size reduction: optional, since size doesn't block PyPI. The experiment
  succeeded (Baseline: 117 MB → 67 MB, 36 MB → 21 MB compressed). To adopt it, add the
  post-build pruning step and the BLAS stub to the build script (B2), check it on macOS
  and Windows (library names differ), and run the tests on the slimmed bundle in CI.
  Less hacky alternatives to look into: whether JuliaC can exclude jll libraries, or
  dropping LinearAlgebra from the dependency tree (hard, because StaticArrays needs it).
- [x] **B7** Threads (2026-10-30, D18). The executable is built with all threads
  (`threads=auto,0`), and the number that work is the number of batches of the parallel
  loops (`_nbatches()`, `src/threads.jl`), also passed to CellListMap (`nbatches`).
  It is set by `-t N`, then the input keyword `nthreads`, then `JULIA_NUM_THREADS`,
  and otherwise all threads are used; `packmol(sys; nthreads)` in Julia. The CLI sets
  BLAS to one thread (Packmol makes no BLAS calls).
  - Timings (s, 32-core machine; the trajectory changes with the number of batches):
    bilayer: 1 core 8.9, `-t 1` 10.0, `-t 4` 4.5, `-t 14` 2.8, all 3.9;
    spherical: 1 core 12.9, `-t 1` 14.0, `-t 4` 8.1, `-t 14` 5.8, all 3.2.
    Single-thread Fortran: 26.3 and 36.7.
  - Pitfall found: with `threads=auto` (one interactive thread by default in Julia ≥
    1.12), a run on one core was up to 3 times slower (bilayer 26.3 s instead of 9.0),
    because the idle interactive thread competes for the core. Hence `auto,0`.
  - `auto` respects the CPU affinity (taskset, Slurm cpusets). Unverified: container
    CPU quotas (`docker --cpus`), where `-t` may be needed.
  - Idle threads cost about 10% (`-t 1` with 32 threads against one thread on one core).
---

## Workstream C: distribution

- [ ] **C1** Release archives on GitHub Releases: `packmol-<ver>-<os>-<arch>.tar.gz`
  (or `.zip` on Windows) with SHA256 sums, produced by B3. Package managers consume these
  archives. Also a rolling `dev` pre-release, replaced on every push to `main`. Link both
  from the installation page of the docs (D20); the `releases/latest/download/<name>`
  URLs need asset names without the version (or a stable alias next to the versioned one).
- [ ] **C2** Python wheels (`pip install packmol`, `uvx packmol`). Reuse the Fortran
  repo's layout: a hatchling package whose `cli.py` execs the bundled binary. The wheel
  carries the prebuilt bundle (no compilation in the build hook). Platform tags:
  `manylinux_2_28_x86_64`/`aarch64`, `macosx_*_arm64`/`x86_64`, `win_amd64`. Size is fine
  (about 36 MB per wheel, or 21 MB slimmed). Wheels can't hold symlinks, so the bundle's
  `lib/` must keep only the soname files the loader looks for (`libjulia.so.1.13`, ...),
  copied as real files.
  The runtime libraries can't be separate PyPI dependencies: they are the exact build of
  one Julia version, and PyPI has no shared package for them. A `packmol-runtime` wheel
  would only split the download without saving anything, so bundling is the way.
- [x] **C3** The Python packaging lives in this repository (D20). The PyPI project
  `packmol` already exists (Fortran wheels, latest 21.2.3, published from m3g/packmol):
  this repository must be added as a trusted publisher of that project, and the new
  versions must sort after 21.2.3 (C5).
- [N] **C4** conda-forge, Homebrew, AUR, distribution packages: left to third parties
  (D20). They are possible by repackaging the C1 archives.
- [ ] **C5** Versioning scheme: binary version, Packmol.jl version and PyPI version.
  Reported by `packmol --version` (A1.5).

---

## Workstream D: features of Packmol.jl that are not in the executable

Document these on a docs page ("Packmol executable vs. Packmol.jl"). Initial list:

- [ ] **D-doc** Write the page. Candidate contents:
  - **Recipes** (`SolutionBoxUS`, `SolutionBoxUSC`, `SolutionBoxUWI`, `Membrane`, with
    Unitful concentrations and `DensityTable`) are a Julia API only. *Option:* expose them
    as CLI subcommands later (`packmol solvate ...`), which would need a trim-safe argument
    parser and maybe no Unitful in that code path.
  - **Julia construction API** (`structure_type(...)`, `PackmolSystem(...)`): no input file needed.
  - **User-defined constraints** (new subtypes of `Constraint`). A compiled binary can't be
    extended; it supports only the constraints built into it, though all of them, including
    gaussian/wave/exponential, are available through the input file.
  - **In-memory interop**: `get_atoms(::PackmolSystem)` with PDBTools / MolSimToolkit /
    ComplexMixtures, and restarting from `Atom` vectors.
  - **`run_packmol_legacy`** (Packmol_jll): Julia-only, and to be dropped.
  - **Interactive file dialog**: the binary has it (`pick_file`), through OS tools.
  - **Differences from Fortran Packmol** (or a separate migration page): restart files
    are coordinate files (D13), `constrain_rotation` constrains the angle between an axis
    and a vector of the molecule (D15), `centerofmass` uses atomic masses (D16),
    `pbc a b c` (D11), the distance criterion (D12), deprecated keywords (D4).
  - Check whether **D=2 packing** and **other number types** (`Float32`) are reachable
    from the CLI (today it is hard-wired to `Val(3), Float64`).

---

## Session log

- **2026-10-29** (Claude): Audit and baseline. Compared the Fortran `getinp.f90` keywords
  with the Julia parser, parsed all Fortran test inputs and examples, built the trimmed
  binary, and ran it against Fortran on 8 tests. Found the CLI-contract gaps (stdin, `-o`,
  exit codes), XYZ support, the unrecognized legacy keywords, inline comments and `nloop`.
  Wrote this plan. Nothing in `src/` changed.
- **2026-10-29** (Claude, same session, after the maintainer's decisions D2/D4–D7):
  - Done: A2.5 (AMBER TER), A3.1 (inline comments), A3.2 (`nloop`/`maxit` from the
    input), A3.3 (deprecation warnings), with tests in `test/data_structures/PackmolSystem.jl`
    and `test/io/write_output.jl`, and updated docs.
  - B6 experiment: slimming the bundle works (Baseline).
  - Then, also done: A1.1–A1.5 (stdin, `-o`, exit codes, `_FORCED`, `-h`/`--version`),
    `segid`, `movefrac`, `movebadrandom`. `read_packmol_input(::IO)` was added. Tests
    are in `test/CLI.jl`.
  - Trim gotchas found along the way:
    - `readlines` is not trim-safe (untyped `ondone` closure), so the parser uses a
      `readline` loop (`_read_lines`).
    - The exception in a `catch` block is `Any`, so pass it to helpers marked
      `@nospecialize`, and narrow `e.msg::AbstractString` to `String` before printing.
    - `Base.Filesystem.File(RawFD(0))` fails on pipes (`lseek`), so stdin is read with
      libc `read`.
  - Windows paths (`_isatty`/`_read` from msvcrt) are untested (B4).
  - All of this went into one PR.
- **2026-10-29** (Claude, same session): XYZ (A2.1) and CONECT (A2.4), using the
  PDBToolsBase `connectivity-title-xyz` branch (1.1.0-DEV, not yet released). Fixed a bug
  where a fixed structure without `number` had 0 molecules. Another trim gotcha: build
  `Structure` metadata with the positional `MetaData(unitcell, connectivity, title)`,
  because a keyword call with the `Union`-typed connectivity is unresolved.
- **2026-10-29** (Claude, same session):
  - Output is always written with `write_pdb`, and TER records are added in a second
    pass.
  - Fortran's benzene2 test showed that PDBToolsBase dropped repeated CONECT bonds
    (bond orders). Fixed in m3g/PDBToolsBase.jl#7, folded into the pending 1.1.0
    registration.
  - Test files moved to MolSimTestFiles (`packmol-v1`, lazy artifact, `TestFiles`
    test module).
  - The binary was run on all inputs of the set (see A4.1 for the 3 gaps).
  - Then the 3 gaps were fixed:
    - A fixed molecule outside the cell is an input error, which also exposed the
      `_collect_fixed_positions` bug.
    - `pbc a b c` now uses Fortran's convention (D11).
  - The Fortran test suite is a Julia test. Test files are at `packmol-v2`.
- **2026-10-29** (Claude, same session): decisions D13–D16 after PR #43 (fixed molecules,
  Fortran restart angles). Removed the raw restart format and the angle conversion of
  PR #43: `restart_to` writes coordinates (PDB/XYZ) and `restart_from` reads any
  coordinate file (`_check_restart_files` at parse time and in `packmol()`). Free
  molecules are now aligned to their principal axes (`set_reference_coordinates!`,
  dimension-generic, proper and idempotent, called from `_center_reference_coordinates!`):
  it existed but was never called. A5.5 and A5.9 closed.
- **2026-10-29** (Claude, same session): the maintainer revised D15 (`constrain_rotation`
  as the angle between an axis and a vector between two atom sets; A5.11, syntax to
  decide). Seed study of `mixture_pbc` (A5.10). The test `structure radius larger than
  tolerance/2` exposed that convergence ignored the radii: D12 refined to a per-pair
  overlap (`InteratomicDistanceFG.max_overlap`).
- **2026-10-29** (Claude, same session): investigated `mixture_pbc` (the maintainer: any
  input taking more than a few seconds signals a bug). Variants showed only PBC + planes
  was slow; the cause was the wrap of centers of mass on the upper faces of the cell
  (A5.10, fixed). The PBC confinement must stay (maintainer): it gives non-periodic
  constraints their meaning. Regression tests in
  `test/function_and_gradient/interatomic_distance_fg.jl`; the seed pin was removed.
- **2026-10-29** (Claude, same session): PR #43 merged. New branch `rotation-constraints`:
  the `constrain_rotation` of D15 (A5.11), with the block syntax approved by the
  maintainer; tests for the penalty, its gradient (finite differences on the whole
  objective), the parser and packings measured from the output. Plan revised (Status
  section, stale items). Then the large test (A4.3).
- **2026-10-29** (Claude, same session): **large test** (A4.3), commit d7c6d81 (branch
  `rotation-constraints`), Julia with 32 threads (`-t auto`), not the binary: converged
  ("Success!") at loop 35; packing 1977.6 s (33.0 min), wall clock 33:09 including
  startup and compilation, peak memory 6.1 GB, output PDB 209 MB.
- **2026-10-30** (Claude, same session): merged PR #44. `xygauss`, `writeout` and
  `writebad` not needed (D17). Branch `build-ci`: `juliac/` build script and executable
  tests (B2, A4.1: 81 checks pass locally in 15 s, build 16 s), `BuildExecutable.yml`
  (B3) on five platforms.
- **2026-10-30** (Claude, same session): PR #45 CI green on all five platforms (B3,
  B4 findings). Found that the executable is single-threaded by default (B7, decision
  needed). Local timings in this log were run with `JULIA_NUM_THREADS=14` set in the
  maintainer's shell.
- **2026-10-30** (Claude, same session): B7 on branch `threads` (D18): `-t`, `nthreads`,
  `_nbatches()` for every parallel loop and CellListMap, BLAS one thread, build with
  `threads=auto,0` (the interactive thread made single-core runs 3 times slower).
- **2026-10-30** (Claude, same session): merged PRs #45 and #46. Branch `tutorials`:
  tutorial pages for the five examples of the Packmol site (command line and Julia),
  `atom_constraints` keyword of `structure_type` (the Julia API had no per-atom
  constraints), D19 (box/cube positional constructors), A5.12 (cube corner fix).
- **2026-10-30** (Claude, same session): merged PRs #47–#50 (tutorials, docs layout,
  MolSimToolkit port). D20 (distribution decisions); B5 finding (glibc 2.34).
- **2026-10-30** (Claude, same session): bug fixes found while writing tutorials, merged:
  `outside box/cube` penalized points outside in only one dimension (#54); `packmol` stops
  when the initial approximation does not satisfy the constraints, unless
  `continue_failed_initial_packing` (#55, exit code 171); fixed atoms on the faces of the
  periodic cell rejected by the rounding of PDB coordinates (#57). Tutorials: new VMD
  images and `check`-first advice for the classic examples (#56); protein in
  glycerol/water, water drop on graphite, gases in chabazite (triclinic) (#58). Scripts in
  `dev/tutorials/` (`render.jl` renders with VMD offscreen).
  **In progress (uncommitted, in the working tree):** the POPC membrane in water/ethanol
  (`docs/src/tutorials/files/membrane/`; the extended POPC conformation is residue 10 of
  ComplexMixturesExamples example 3) does not pack: lipids alone fail at 64–75 Å²/lipid,
  even with head/tail planes 47% of the lipid length apart; Fortran does worse; our
  `Membrane` recipe fails with this POPC too. To study (less dense membranes; chains as
  separate molecules; flexibility). The vesicle and lipid nanotube inputs are written but
  not run; polymer in DMF has files only. Not started: the Zika virion tutorial (its input
  packs NaW twice; the last block, under "Chloride", should probably be ClW), the
  50%-rule SVG and recipe default, and the very large (10⁶–10⁷ atoms) system experiment.
- **2026-10-30** (Claude, same session): bug fixes found while writing tutorials, merged:
  `outside box/cube` penalized points outside in only one dimension (#54); `packmol` stops
  when the initial approximation does not satisfy the constraints, unless
  `continue_failed_initial_packing` (#55, exit code 171); fixed atoms on the faces of the
  periodic cell rejected by the rounding of PDB coordinates (#57). Tutorials: VMD images
  and `check`-first advice for the classic examples (#56); protein in glycerol/water,
  water drop on graphite, gases in chabazite (triclinic) (#58); Zika virion (files in the
  `packmol_zika-v1` release of MolSimTestFiles), zeolite with a live PDBTools view. The
  POPC membrane tutorial does not pack yet (see the memory notes / tutorials-status):
  lipids alone fail at 64–75 Å²/lipid with an extended POPC and permissive planes.
