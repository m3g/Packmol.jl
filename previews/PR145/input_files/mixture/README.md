# packmol_tutorials

The structure files of the tutorials of [Packmol.jl](https://github.com/m3g/Packmol.jl)
(`docs/src/tutorials`), in a single directory: each file is used by every tutorial that
needs that molecule (one water, one POPC, one gold core, ...). The input files and the
Julia scripts of the tutorials are in the documentation of Packmol.jl; the larger
tutorials have their own sets (`packmol_zika`, `packmol_exosome`, `packmol_lnp`, and the
Bentopy models of the minimal cell).

- `water.pdb`: TIP3P water (all-atom tutorials); `water.xyz`, `chlor.xyz`, `t3.xyz`: the
  water–chloroform interface (xyz input); `wt4.pdb`: the WT4 water of SIRAH (condensate).
- `popc.pdb`: POPC, as in the `packmol_lipids` set.
- `thyroid_receptor.pdb`: the ligand-binding domain of the thyroid hormone receptor
  (Martínez et al., J. Med. Chem. 49, 23, 2006), solvated protein tutorial.
- `subtilisin.pdb`: subtilisin, protein in glycerol and water tutorial.
