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
  Here a block without `segid` keeps the segment names of its own PDB file.
- **D10 (2026-10-29): `movebadrandom` is deprecated** (warned and ignored). Julia's `movebad!` always
  places bad molecules at random positions (within the per-type center-of-mass bounds,
  with a short constraint fit), which is Fortran's `movebadrandom` behavior. Fortran's
  default instead places them next to a randomly chosen good molecule. That mode is not
  implemented; add it only if benchmarks show it helps.

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
- [ ] **A2.1** `filetype xyz`: XYZ structure input and XYZ output (used in
  `examples/interface.inp`). It is implemented in PDBToolsBase first (D5), then wired here.
- [?] **A2.2** `filetype tinker` (and `moldy`): marked "won't be supported" in CLAUDE.md.
  Today `filetype` is parsed and silently ignored, so a Tinker input would be read as PDB.
  At least fail with a clear error for any `filetype` other than pdb (and xyz once A2.1
  lands).
- [ ] **A2.3** `writecrd` (CHARMM CRD output, Fortran-only today). Decide: implement or reject clearly.
- [ ] **A2.4** CONECT (`connect`, `ignore_conect`, `non_standard_conect`, `hexadecimal_indices`):
  needs PDBToolsBase support first (D6). The three option keywords are warned and ignored
  for now.
- [x] **A2.5** AMBER TER (`add_amber_ter`: a TER after every molecule; `amber_ter_preserve`:
  the template's TER records are repeated in every copy), as in Fortran `output.f90`.
  Done 2026-10-29 (`_write_pdb_with_ter` in `src/io/write_output.jl`).

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
  `precision` and `iprint*`. `segid` (per structure) and the CONECT options are
  "not yet implemented" warnings.
- [x] **A3.3b** `segid` (D9).
- [ ] **A3.4** `over xygauss` / `below xygauss` (Fortran's Gaussian surface): map onto
  the `AboveGaussian`/`BelowGaussian` constraints (check that the parameterization is
  equivalent).
- [x] **A3.5** `movefrac` is read into `sys.movefrac` (the `packmol()` keyword argument
  takes precedence). `movebadrandom` is deprecated (D10).

### A4. Validation
- [ ] **A4.1** Port Fortran's `testing/` suite (runtests.jl validates outputs with
  CellListMap; test.sh, test_cli.sh, test_failed.sh) to run **against the compiled
  binary**, with the inputs adapted for D1. This is the drop-in acceptance test.
- [ ] **A4.2** Benchmark the binary against Fortran on the larger examples (solvprotein is
  slower single-threaded), and record the results here.

---

## Workstream B: standalone executable build (and CI)

- [x] **B1** Trim-safe build works locally on Linux (see Baseline).
- [ ] **B2** Build script in the repo (`build/` or `juliac/`: `app.jl` + a script running
  juliac), so that CI and humans run the same thing. Pin the Julia version used to build.
- [ ] **B3** CI workflow `BuildExecutable.yml`: a matrix with ubuntu-latest,
  ubuntu-24.04-arm, macos-latest (arm64), macos-15-intel and windows-latest, each building
  the bundle, running smoke tests (A4.1 subset) on the bundle, and uploading it as an
  artifact. Trigger on tags/releases and `workflow_dispatch`, plus a cheap Linux run on PRs
  to catch trim regressions (any new non-trim-safe call breaks the build).
- [ ] **B4** **Windows and macOS spikes.** It is unverified whether `juliac --trim` works
  on them with 1.13. On Windows, look at the toolchain JuliaC needs, and at `pick_file`'s
  PowerShell path. On macOS, look at code signing/notarization (an unsigned binary is
  quarantined when downloaded through a browser, though pip/conda installs are fine).
  Record the findings here.
- [ ] **B5** Linux portability: the build host's glibc sets the minimum (≥ 2.34 today).
  Build in a manylinux_2_28 container (or older) for wider coverage, matching the wheel tag.
- [~] **B6** Size reduction: optional, since size doesn't block PyPI. The experiment
  succeeded (Baseline: 117 MB → 67 MB, 36 MB → 21 MB compressed). To adopt it, add the
  post-build pruning step and the BLAS stub to the build script (B2), check it on macOS
  and Windows (library names differ), and run the tests on the slimmed bundle in CI.
  Less hacky alternatives to look into: whether JuliaC can exclude jll libraries, or
  dropping LinearAlgebra from the dependency tree (hard, because StaticArrays needs it).
- [ ] **B7** Threads: document `JULIA_NUM_THREADS` and/or add a `-t N` flag. Consider
  defaulting to several threads.

---

## Workstream C: distribution

- [ ] **C1** Release archives on GitHub Releases: `packmol-<ver>-<os>-<arch>.tar.gz`
  (or `.zip` on Windows) with SHA256 sums, produced by B3. Package managers consume these
  archives.
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
- [?] **C3** Where the Python packaging lives: in this repo, or in m3g/packmol (which then
  becomes the "distribution" repo, downloading the binaries from Packmol.jl releases).
  This needs a decision (see D2).
- [ ] **C4** conda-forge recipe (repackaging the C1 binaries, or building from source with
  juliac, which is harder), then Homebrew, AUR, etc. Lower priority, after C1/C2.
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
