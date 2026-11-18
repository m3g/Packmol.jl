# Structure blocks

```
structure <pdb file>
    number <n>
    [fixed x y z beta gamma theta]
    [connect no]
    [center | centerofmass]
    [radius <value>]
    [radscale <value>]
    [resnumbers <0|1|2|3>]
    [chain <c> | changechains]
    [segid <id>]
    [constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>
      [atoms <i> ...
       to <j> ...
     end constrain_rotation]]
    [restart_from <file>]
    <constraint lines...>
    [atoms <indices...>
        <per-atom overrides...>
    end]
end structure
```

| Keyword | Meaning |
|:--|:--|
| `number <n>` | Number of copies of this molecule to pack (required). |
| `fixed x y z beta gamma theta` | Don't pack this molecule — place it once at this position and (Euler angle) orientation and leave it fixed. `number` is optional for a fixed structure (it is a single molecule). See [Fixed molecules](../constraints/fixed.md). |
| `connect no` | Don't write the connectivity (`CONECT` records) of this structure to the output. |
| `center` / `centerofmass` | Only valid together with `fixed`: re-center the molecule (geometric center or mass-weighted center of mass) before applying the fixed position/orientation. Without them, the molecule is rotated about the origin of its coordinates, and then translated. (Fortran Packmol uses unit masses for `centerofmass`, so its center of mass is the geometric center.) |
| `radius <value>` | Atom radius for every atom of this structure, overriding the default of `tolerance/2`. |
| `radscale <value>` | Scale of the radii of the atoms of this structure at the start of the packing, decaying to 1 in `radscale_decay_loops` loops (see [Global keywords](@ref)). Overrides the global `radscale` and the default: `5.0` for free molecules of more than 20 atoms, `1.0` for the other free molecules and for restarted ones (`restart_from`), and the largest of those for fixed molecules. |
| `resnumbers <n>` | Residue numbering scheme for the output — see [Residue numbering](@ref) below. |
| `chain <c>` | Fix every molecule of this type to chain identifier `<c>` in the output, instead of the default auto-assigned one. |
| `changechains` | Alternate between two auto-assigned chain identifiers every other molecule, instead of one chain for the whole type. Not compatible with `chain`. |
| `segid <id>` | Segment identifier (at most 4 characters) for every atom of this type, in columns 73-76 of the output PDB file. Without `segid`, the segment identifiers of the structure's own PDB file are kept (blank if it has none). |
| `constrain_rotation <axis> <center_deg> <halfwidth_deg>` | Keep the angle between `<axis>` (`x`, `y` or `z`) and a vector of the molecule within `<center_deg> ± <halfwidth_deg>` degrees: the longest axis of the molecule, or, followed by `atoms <i> ...`, `to <j> ...` and `end constrain_rotation` lines, the vector between the centers of two sets of atoms. See [Rotation constraint](@ref). |
| `restart_from <file>` | Restart only this structure type's molecules — see [Restarting a run](@ref). |
| `atoms <indices> ... end` | Restrict the enclosed lines (`radius`, constraints) to the listed 1-based atom indices instead of the whole molecule. |

Constraint lines (`inside box`, `outside sphere`, `above plane`, ...) are
documented on the [Constraints](../constraints/index.md) page — any number of them can
appear in a structure block, optionally inside an `atoms ... end` block.

## Several conformations of a molecule

A structure file with several models (`MODEL` ... `ENDMDL` records of a PDB file)
provides several conformations of the same molecule: the `number` of molecules is
divided among them as evenly as possible, and every constraint of the structure block
applies to all of them. The models must have the same atoms, in the same order (the
atom indices of `atoms` blocks refer to all of them).

```
structure popc_conformations.pdb # 30 models
  number 60                      # 2 molecules of each conformation
  atoms 1
    above plane 0. 0. 1. 19.5
  end atoms
  atoms 88 131
    below plane 0. 0. 1. 6.
  end atoms
end structure
```

Molecules are rigid in Packmol: several conformations of a flexible molecule pack
more densely than copies of a single one, since the molecules fill the space between
each other better. Lipid membranes, for example, pack at their actual density from
conformations taken from an equilibrated membrane, but not from copies of a single
extended lipid. In the output, the molecules of all conformations are numbered as
those of a single structure type. A fixed molecule has a single conformation, and
`restart_from` is not available inside a structure block with several conformations
(use it for the whole system).

## Polymers

A structure block with `repeat` blocks is a polymer (Packmol.jl only): the structure file
holds the whole polymer (one model per conformation), and the `repeat` blocks give its
sequence of monomers, each with its number of atoms, in the order of the atoms of the
file. Each monomer is packed as a rigid molecule, and consecutive monomers are held
together by **links**: harmonic restraints between an atom of a monomer and an atom of
the next one, at their distance in the structure. The chains are thus flexible: the
monomers turn about the junctions, while the bonds and angles of the junctions are kept
by the links.

```
structure peptides.pdb      # the polymer: 10 models, 10 conformations
  number 123                # polymers
  inside sphere 0. 0. 0. 100.
  repeat 20                 # 20 times the monomers of the block
    monomer 18              # 18 atoms
      link 17 1             # atom 17 of this monomer to atom 1 of the next
      link 18 1
      link 18 2
  end repeat
end structure
```

| Keyword | Meaning |
|---|---|
| `repeat <n> ... end repeat` | `n` repeats of the monomers of the block. A sequence can have several blocks, one after the other (block copolymers, end groups with `repeat 1`, or an irregular sequence as a single block with `repeat 1`). |
| `monomer <natoms>` | A monomer of the sequence, of `natoms` atoms (the next `natoms` atoms of the structure file). |
| `link <a> <b> [d0 [k [w]]]` | A link from atom `a` of the monomer (the last `monomer` line) to atom `b` of the next monomer (atom indices local to each monomer). Its penalty is `k (\|d - d0\| - w)²` outside `d0 ± w`: `d0` (Å) is by default the distance of the two atoms in each conformation, `k` 10, and `w` 0 (no flat bottom). The links of the last monomer of a polymer are ignored, unless it is `cyclic`. |
| `cyclic` | The last monomer is linked to the first one, with the links of the last monomer. |

Two links per junction fix the bond and the angles of the junction, and leave free the
rotation about the bond: as many links as needed to keep the geometry of the junction
(three, here, from the last two atoms of a monomer to the first two of the next one).

- **The structure of the polymer.** The monomers have the geometry of each conformation
  of the file: the number of polymers is divided among the conformations, as for other
  molecules (see above), and each polymer is built from the monomers of its own.
- **Packing.** The polymers are first packed as rigid molecules, in their conformations,
  and then their monomers are released and the packing goes on with the flexible chains.
  (Chains released overlapping each other end up threaded through each other's monomers,
  and do not separate.) Polymers too dense to pack as rigid molecules are packed at a
  lower density, and compressed (see [Packing by compression](@ref)).
- **Close atoms.** The pairs of atoms of different monomers (up to three monomers apart)
  closer than the tolerance in a conformation (the bonded atoms of a junction, and their
  neighbors) are left out of the distances of the polymers of that conformation.
- **The convergence criterion** includes the links: every link within `d0 ± w` to the
  precision of the tolerance (`tolerance_precision`).
- **Defaults.** The radii of the atoms of polymers are not inflated at the start of the
  packing (`radscale` 1.0: inflated atoms push consecutive monomers apart, against the
  links); movebad moves the monomers of a bad polymer by small random displacements and
  rotations (a polymer moved as a whole to a random place, in a dense system, lands across
  the others).
- **Output.** Each polymer is written whole, with the atoms of the structure file, in
  their order (with the residue numbers of the file), a chain identifier per polymer (or
  the `chain` of the block), and the bonds of the file (`CONECT` records); with periodic
  boundary conditions, a polymer is not wrapped across the boundaries of the cell (it is
  whole, its first monomer in the cell).
- **Not supported**: fixed polymers, `restart_from`, `constrain_rotation` and
  `voxel_place`.

Polymers are a case of flexible molecules, as molecules with [Rotatable bonds](@ref):
the monomers are the rigid fragments, given by their numbers of atoms, and the links are
written explicitly (the structure files of coarse-grained polymers rarely have bonds).

In Julia, the polymers are given by [`polymer_structure`](@ref), with the sequence made of
[`monomer`](@ref)s. The tutorial [Biomolecular condensate](@ref) packs a condensate of
polypeptides, in a homogeneous mixture and as a droplet.

## Rotatable bonds

Any molecule can be made flexible about some of its bonds (Packmol.jl only): with
`rotatable` lines, the molecule is cut at those bonds into rigid fragments (the
connected parts of the bonds that are left), packed as rigid molecules held together by
links. The links of a rotatable bond i–j keep its length and its angles: the bond i–j,
and the bonds of i and of j to their other neighbors in their fragments. The fragments
turn about the bond, and nothing else changes.

```
structure popc.pdb
  number 60
  bond_distance 1.8                                # the bonds, from the distances
  rotatable 25 28  33 45  57 60  42 92  101 104    # pairs of bonded atoms
  atoms 1
    above plane 0. 0. 1. 28.
  end atoms
end structure
```

| Keyword | Meaning |
|---|---|
| `rotatable <i> <j> [<i> <j> ...]` | Rotatable bonds: pairs of bonded atoms (indices of the structure file), on one or several lines. A bond in a ring cannot rotate (an error). |
| `bond_distance <d>` | The bonds of the molecule are the pairs of atoms closer than `d` Å in the structure (except pairs of hydrogens), instead of the `CONECT` records of the file. `1.8` finds the bonds of organic molecules with hydrogens (the hydrogens of a CH₂ are 1.75 Å apart, as close as some bonds). |

The rotatable bonds are chosen by hand: each one adds a rigid fragment, with six
variables. Choose them where the molecule must bend to fit: a chain that has to fold
needs to turn about many of its bonds (see the example below). Without
`bond_distance`, the bonds are those of the `CONECT` records of the structure file (an
error if it has none). The pairs of atoms of different fragments up to three bonds apart
are left out of the distances, and the rest is as for [Polymers](@ref): the molecules
are first packed as rigid molecules, then the fragments are released; each molecule is
written whole, in the order of the atoms of the file; the convergence criterion includes
the links.

For example, the 120 POPC lipids of the [Lipid bilayer](@ref) tutorial (64 Å² per lipid)
in thinner membranes: the head (atom 1) within a band of 3 Å above `H` (and below `-H`),
the ends of the chains between -2 and 6 Å (and 2 and -6 Å), and the water beyond the
heads, at the density of the tutorial. The lipid of the structure is extended, 28.9 Å from
the head to the ends of the chains:

| Heads at (Å) | Rigid | 5 rotatable bonds | 28 rotatable bonds | 28, with `compress 0.5` |
|---|---|---|---|---|
| 24 to 27 | packs | packs | | |
| 20 to 23 | does not pack (f = 65) | does not pack (f = 8) | packs (36 s) | packs (76 s) |
| 17 to 20 | does not pack (f = 384) | does not pack (f = 36) | does not pack (f = 21) | 96% of the density |

The 28 rotatable bonds are all the single C–C bonds of the chains (but those of the
methyl groups) and the C1–C2 bond of the glycerol:

```
  bond_distance 1.8
  rotatable 25 28  33 45  45 48  48 51  51 54  54 57  57 60  60 63  65 67  67 70  70 73  73 76  76 79  79 82  82 85
  rotatable 42 92  92 95  95 98  98 101  101 104  104 107  107 110  110 113  113 116  116 119  119 122  122 125  125 128
```

In the packed membrane with the heads at 20 to 23 Å, the torsions turned by up to 81°,
and the heads are 23 to 29 Å from the ends of the chains. With `compress 0.5` (see
[Packing by compression](@ref)), the membrane is first packed 26% thicker, where the
rigid lipids fit, and then thinned with the lipids flexible: the rigid lipids, compressed
alone, stop at 82% of the density.

In Julia, the molecules with rotatable bonds are given by [`flexible_structure`](@ref).
