# Contributing to Packmol.jl

Contributions are welcome: bug reports, examples that fail to pack, documentation,
new constraints, recipes and improvements of the packing. These guidelines apply to
everyone, the maintainers included.

## Reporting problems

Open an [issue](https://github.com/m3g/Packmol.jl/issues) with the input file and the
structure files needed to run it (or the smallest system that shows the problem), the
version of Packmol.jl (`packmol --version`, or `pkg> st Packmol`), and the output. A
system that does not pack is a useful report: say how long it ran and what the last loops
of the output show.

## Setting up

```bash
git clone https://github.com/m3g/Packmol.jl
cd Packmol.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Julia 1.13 or later is required. The layout of the code is described in `src/README.md`,
and each directory of `src/` has a README describing its files.

## Pull requests

- One topic per pull request, from a branch (not from `main`).
- Every pull request that changes the package gets an entry in `CHANGELOG.md`, in the
  section of the version under development (`Version X.Y.Z-DEV`), with a badge:
  `FEATURE`, `ENHANCEMENT`, `BUGFIX`, `DEPRECATION`, `BREAKING` or `EXPERIMENTAL`.
  Changes of the documentation only get an `INFO` entry. The `changelog` workflow checks
  it; changes that need no entry (e.g. of the workflows) get the `Skip-Changelog` label.
- All the checks of the pull request must pass (see below) before it is merged.

## The tests that must pass

### 1. The test suite

```bash
julia --project -t auto -e 'using Pkg; Pkg.test()'
```

The tests are `@testitem`s (TestItems.jl) in `test/`, whose layout mirrors that of
`src/`: the tests of `src/constraints/boxes.jl` are in `test/constraints/boxes.jl`. New
code gets tests there. The test data files (structures and input files) are not in this
repository: they are the `packmol` set of
[MolSimTestFiles](https://github.com/m3g/MolSimTestFiles), a lazy artifact, accessed in
the tests with `setup=[TestFiles]` and `TestFiles.test_file(name)`. To add test files,
a new release of that set is published, and `Artifacts.toml` updated.

The CI runs the suite on Linux, macOS and Windows, with Julia 1.13 and the pre-release
(`CI.yml`), and with the lowest versions of the dependencies allowed by `Project.toml`
(`Downgrade.yml`). Changes of the documentation and Markdown files only skip it.

### 2. The executable (`juliac --trim`)

The `packmol` command line program is compiled with `juliac --trim=safe`, which requires
every call reachable from it to be resolved at compile time. Changes of the engine (the
code reachable from `src/CLI.jl`) must keep it compiling:

```bash
julia --project=juliac juliac/build.jl            # the bundle, in juliac/build
julia --project=juliac juliac/test_executable.jl  # runs the examples with it
```

The CI builds and tests it on five platforms (`BuildExecutable.yml`). The rules that keep
the code trim-compatible (printing with `_print`/`_println`/`@_printf`, `@_sync`, no
closures over reassigned variables, ...) are listed in the section "`juliac --trim`
compatibility" of `CLAUDE.md`.

### 3. The documentation

```bash
julia --project=docs --threads=auto docs/make.jl
```

The documentation runs its examples, and packs the systems of the tutorials shown live, so
it is also a test (`Documentation.yml`). The build fails if an example fails, a
reference is broken, or a docstring of the package is not in the pages (new functions
get their docstring in the reference of the Julia API).

### 4. The large tests, for changes of the packing heuristics

Changes of the packing itself (movebad, the initial placement, the compression, the
optimizer, the radii, the active set) can help some systems and hurt others, and the test
suite packs only small systems. Such changes are validated on the large examples (the
tutorials, from seconds to minutes each), against `main`, with several seeds, as packings
with several threads are not deterministic:

```bash
git checkout main
julia --project=dev/large_tests -t 14 dev/large_tests/run.jl --seeds 1:3 --skip exosome --label main --results validation.tsv
git checkout <branch>
julia --project=dev/large_tests -t 14 dev/large_tests/run.jl --seeds 1:3 --skip exosome --label new --results validation.tsv
julia dev/large_tests/compare.jl validation.tsv main new
```

`compare.jl` flags the cases that converged less often with the change, or got slower by
more than 30%. A change is accepted if no case gets worse: report the comparison in the
pull request. See `dev/large_tests/README.md`. Run the cases one at a time, with nothing
else heavy on the machine, or the timings are meaningless. Some changes that looked
good on one system were rejected this way: local moves in movebad, for instance, helped no
case and broke the polymer-capped gold nanoparticle.

New tutorials with non-trivial packings get a case in `dev/large_tests/run.jl`.

## Code

- Dimensions are parameters (`D`, the 3D case being the target), coordinates in Å.
- Errors are `ArgumentError`s with the message padded by blank lines, every line written
  literally in the string:

  ```julia
  throw(ArgumentError("""\n
      use_compression is not supported with restart_from.

  """))
  ```

- Comments explain why the code is as it is, in particular the alternatives that were
  tried and failed (with the system and the result), so that they are not tried again.
- Input keywords and the Julia API go together: a new input keyword gets its keyword
  argument (in `packmol`, `structure_type` or `PackmolSystem`), and both are documented.

## Documentation

- The documentation is in `docs/src/`, a directory per tab of the top menu, a page per
  topic. New constraints and periodic cells get a page ending with an "Example" section,
  with the structure block of an input file and the same system packed in an `@example`
  block and shown with `PDBTools.visualize`.
- The pages describe Packmol.jl as it is: comparisons with the legacy (Fortran) Packmol
  are kept to `docs/src/migrating.md`.
- Julia code in the documentation is written for readability: one assignment per line,
  no compact patterns that only save lines.
- Tutorials: their input files are files of `docs/src/tutorials/files/<tutorial>/`, shown
  in the pages, not written into them; each starts with sections on the system (its
  structures, regions and numbers of molecules), then sections named after what they
  build, and a last section, "Julia script", with the same system in Julia. Data that
  cannot be redistributed is not added: the page points to its source.

## License

Packmol.jl is distributed under the MIT license (`LICENSE`). Contributions are accepted
under the same license.
