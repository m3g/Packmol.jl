Packmol.jl Changelog
====================

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 100.1.0-DEV
-------------------
- ![FEATURE][badge-feature] Polymers (Packmol.jl only): a structure block with `repeat`
  blocks (`repeat <n>`, `monomer <natoms>`, `link <a> <b> [d0 [k [w]]]`, `cyclic`), or
  `polymer_structure` and `monomer` in Julia, packs flexible chains of rigid monomers held
  together by links between their atoms, or by bonds about which they turn (`rotatable
  <a> <b>` in a monomer, with the bonds from CONECT records or `bond_distance`). The
  polymers are packed first as rigid molecules, then released.
- ![FEATURE][badge-feature] Rotatable bonds (Packmol.jl only): `rotatable <i> <j> ...` in a
  structure block (`flexible_structure` in Julia) cuts the molecule into rigid fragments at
  those bonds, linked such that they turn about the bonds only. The bonds are those of the
  CONECT records, or computed from the distances (`bond_distance <d>`).
- ![FEATURE][badge-feature] `use_compression [<fraction> [<step>]]` (Packmol.jl only; a keyword
  argument in Julia): packs the system at a fraction of its density (default 0.5), with the
  constraints and the periodic cell expanded, and compresses it in adaptive steps to the
  density of the input, for systems too dense to pack from a random placement. A step that
  does not converge is relaxed (expanded below the last converged density, keeping its
  rearrangements, and packed) before the compression resumes with half the step; a packing
  at the low density that does not converge is tried again at half that density, up to
  three times. Fixed molecules are scaled with the geometry (the coordinates of their atoms
  and their radii), and restored at the end. All periodic cells are supported; not with
  `restart_from` or `voxel_place`.
- ![FEATURE][badge-feature] `rotatable all` (Julia: `flexible_structure(...; rotatable = :all)`):
  every rotatable bond of a molecule, found from its bonds: not in rings, between heavy
  atoms with other heavy neighbors, and single for sure from the numbers of bonds of the
  atoms (bonds between two carbon or nitrogen atoms with fewer than four bonds do not
  turn, but those between aromatic rings). For all-atom structures with all their hydrogens.
- ![INFO][badge-info] `CONTRIBUTING.md`: the guidelines for contributions, and the checks a
  pull request must pass (the test suite, the executable, the documentation, and the
  validation on the large tests of changes of the packing heuristics).
- ![INFO][badge-info] Documentation: the tutorial "Lignin particle in acetone and water": 70
  flexible oligomers of spruce lignin (built with LigninBuilder), with `rotatable all`,
  packed by compression at 1.30 g/cm³, then solvated; and dispersed in the solvent.
- ![INFO][badge-info] Documentation: the flexible molecules (rotatable bonds, polymers) and the
  packing by compression in "New in Packmol.jl" of the home page.
- ![INFO][badge-info] Documentation: the tutorial "Polymer-capped gold nanoparticle": 57
  flexible chains of thiol-terminated polystyrene (all atoms) bound to a fixed gold core, in
  toluene, packed by compression around the core.
- ![INFO][badge-info] Documentation: the tutorial "Biomolecular condensate", with chains of
  (VPGVG)₂₀ (SIRAH) in water, as a homogeneous mixture and as a droplet at 400 mg/mL, and the
  pages "Polymers" (Structure blocks) and "Packing by compression".
- ![INFO][badge-info] Documentation: the overview of the tutorials links the classic examples
  to the first paper of Packmol (2003).
- ![INFO][badge-info] Documentation: a "See also" section in the introduction (related packages
  and contact).
- ![ENHANCEMENT][badge-enhancement] The error of a `rotatable` bond whose atoms are not bonded
  says why: not in the CONECT records, or farther than `bond_distance` (with their
  distance). The docs state that the atom indices count from 1, by position in the file.
  `bond_distance` without rotatable bonds is no longer an error (it has no effect).

Version 100.0.0
---------------
The first release of Packmol.jl as Packmol: a new implementation, in Julia, of the Fortran
Packmol (whose last version is 21.2.3), distributed as a Julia package, a standalone
executable, and Python wheels (`pip install packmol`). Its version numbers continue those of
the Fortran Packmol. Changes relative to the Fortran Packmol (see the "Migrating from legacy
Packmol" page of the documentation):

- ![INFO][badge-info] The input files of the Fortran Packmol run unchanged in most cases,
  with `packmol < input.inp`, `packmol input.inp` or `packmol -i input.inp`. The configurations
  differ from those of the Fortran Packmol for the same input and seed (new optimizer, initial
  placement and random number generator).
- ![BREAKING][badge-breaking] The paths of an input file (structures, `output`,
  `restart_from`) are relative to the directory of the input file (with `packmol < input.inp`,
  to the current directory, as before).
- ![BREAKING][badge-breaking] `restart_from` reads coordinate files (PDB, XYZ or mmCIF, for
  example a previous output), not the restart files of the Fortran Packmol; `restart_to` is
  ignored, with a warning.
- ![BREAKING][badge-breaking] `constrain_rotation` constrains the angle between an axis and a
  vector of the molecule (its longest axis, or the vector between two sets of atoms), not the
  Euler angles.
- ![BREAKING][badge-breaking] `centerofmass` uses the masses of the atoms (the Fortran
  Packmol used the geometric center).
- ![BREAKING][badge-breaking] Misspelled keywords in structure blocks are errors.
- ![BREAKING][badge-breaking] The formats `tinker` and `moldy`, and the constraints
  `over xygauss`/`below xygauss` (replaced by `above gaussian`/`below gaussian`), are not
  supported.
- ![BREAKING][badge-breaking] The output that satisfies only the constraints, written when the
  packing does not converge, is `system-FORCED.pdb` for the output `system.pdb` (was
  `system.pdb_FORCED`).
- ![DEPRECATION][badge-deprecation] Keywords of the old algorithm are accepted with a warning
  and ignored: `fscale`, `fbins`, `precision`, `sidemax`, `maxmove`, `nloop0`, `packall`,
  `disable_movebad`, `movebadrandom`, `iprint1`/`2`/`3`, the `short_*` keywords,
  `randominitialpoint`, `optimization_print_level`, `chkgrad`, `writebad`, `writecrd`,
  `non_standard_conect`, `hexadecimal_indices`, and `nloop` in structure blocks. `check` is now
  `check_initial_placement`.
- ![ENHANCEMENT][badge-enhancement] Faster and multithreaded (`packmol -t N`): systems of
  millions of atoms (membranes, vesicles, virions, cells) pack in minutes.
- ![FEATURE][badge-feature] mmCIF files, for structures, outputs and restarts.
- ![FEATURE][badge-feature] Periodic cells: `pbc xmin ymin zmin xmax ymax zmax`, triclinic
  cells (`unitcell a b c α β γ`), rhombic dodecahedra and truncated octahedra
  (`pbc dodecahedral`, `pbc octahedral`).
- ![FEATURE][badge-feature] New constraints: curved surfaces (Gaussian, waves and exponentials,
  and their radial forms), surfaces given by points (`within surface`, `inside surface`), and
  rotation constraints.
- ![FEATURE][badge-feature] `check_initial_placement`, to see where the constraints put the
  molecules before packing; several conformations of a molecule, from files with several
  models; `voxel_place`, the placement of Bentopy, for crowded systems.
- ![FEATURE][badge-feature] Command line: `-o` and `-s` (output and seed, for replicas of a
  system), `--guess-volume` and `--free-volume` (volumes of the system and of the regions of
  its molecules), and a file dialog when run without an input.
- ![FEATURE][badge-feature] Julia API: systems defined without input files, recipes (solvated
  systems, membranes and membrane proteins from target densities and concentrations), custom
  constraints, and the packed atoms in memory (`get_atoms`).
- ![INFO][badge-info] The Fortran Packmol can be run from Julia with `run_packmol_legacy`
  (requires `Packmol_jll`).
