# Large tests

The large examples of Packmol.jl, mostly those of the tutorials, which take from seconds to
hours. They are not part of the test suite (`test/`), and are run on request only, e.g.
before a release. They are plain scripts, not `@testitem`s, so that the test suite does not
pick them up.

From the root of the repository:

```bash
julia --project=dev/large_tests -t 14 dev/large_tests/run.jl --list
julia --project=dev/large_tests -t 14 dev/large_tests/run.jl                 # all cases
julia --project=dev/large_tests -t 14 dev/large_tests/run.jl zika bilayer    # some of them
```

Options: `--seed N` (adds `seed N` to the inputs), `--workdir DIR` (where the inputs are
copied and packed; a temporary directory by default).

The cases are run one after the other (never run two of them, or other heavy jobs, at the
same time: the timings would be meaningless). Each packs a copy of its input directory, and
passes if the packing converged and the atoms of different molecules are at least
`tolerance - tolerance_precision` apart. A line per case is appended to `results.tsv`
(date, commit, Julia version, threads, case, seed, status, loops, function evaluations,
time, minimum distance between molecules, at most the tolerance), which is kept in the
repository to compare versions: commit it after a run on a clean working tree.

The inputs are those of `docs/src/tutorials/files` (cases whose input is missing are
skipped), the Zika virion of the `packmol_zika` set and the exosome of the
`packmol_exosome` set of [MolSimTestFiles](https://github.com/m3g/MolSimTestFiles),
downloaded when first needed (7 MB and 1 MB downloads, `Artifacts.toml`). The exosome
(13.6 million atoms, packed in two steps: the proteins, then the rest) takes about 35
minutes and 24 GB of memory.
