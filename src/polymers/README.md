# polymers

Flexible molecules: polymers (structure blocks with `repeat` blocks, or
`polymer_structure` in Julia) and molecules with rotatable bonds (`rotatable` lines, or
`flexible_structure`), packed as rigid fragments held together by links (harmonic
restraints between atoms of different fragments). Each fragment of each conformation is a
structure type (conformation after conformation), with a molecule for each copy of the
molecule in that conformation. The names keep the words of the polymers (`Polymer`,
monomers): a monomer is a fragment.

| File | Contents |
|---|---|
| `types.jl` | `PolymerLink`, `MonomerSpec`, `PolymerSpec` (the sequence given in the input), `Polymer` (a polymer of a system: its rigid structure types, monomers, links and excluded pairs), `SystemLink` and `PolymerTopology` (the links and excluded pairs of a system), and the helpers to index the monomers. Included before `PackmolSystem`, which holds them. |
| `read_polymer_block.jl` | The `repeat`/`monomer`/`link`/`cyclic` lines of a structure block (`_read_polymer_spec`), and the expansion of a polymer into its monomers (`_add_polymer!`): links resolved (distances from each conformation), excluded pairs (atoms of different monomers, up to 3 apart, closer than the sum of their radii in a conformation). |
| `rotatable.jl` | The molecules with rotatable bonds (`_add_rotatable!`): the bonds (CONECT, or from the distances, without pairs of hydrogens), the fragments (the bond graph cut at the rotatable bonds), the links (each bond and the bonds of its ends to their neighbors), and the excluded pairs (up to three bonds apart); `rotatable all` (`_all_rotatable_bonds`: the bonds certainly single, from the numbers of bonds of the atoms). |
| `topology.jl` | `_set_polymer_topology!` (the links between atoms of the system, the polymer of each monomer molecule, the atoms of excluded pairs), `_excluded_pair` (used by the pair functions), and the indices of the structure types and molecules of the monomers. |
| `link_fg.jl` | The function and gradient of the links, `k (\|d - d0\| - w)²`, for the whole system and for the active subsystem; their largest deviation is part of the convergence criterion. |
| `placement.jl` | `RigidPolymerSystem` (the system with the polymers as rigid molecules), the initial placement (the rigid system placed and packed, then the monomers aligned to their atoms: released), and the movebad of polymers (`_movebad_polymers!`: random displacements of the monomers of bad polymers). |
| `output.jl` | The polymers written whole: wrapped by their centers, in the order of the atoms of their structure, a chain identifier per polymer, and their bonds. |
| `api.jl` | `polymer_structure`, `monomer` and `flexible_structure` (the Julia API). |
