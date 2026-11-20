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

A structure block with `repeat` blocks is a polymer: the structure file
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
| `link <a> <b> [d0 [k [w]]]` | A link from atom `a` of the monomer (the last `monomer` line) to atom `b` of the next monomer (the indices are local to each monomer, from 1: the position of the atom among the atoms of the monomer, in the order of the structure file). Its penalty is `k (\|d - d0\| - w)²` outside `d0 ± w`: `d0` (Å) is by default the distance of the two atoms in each conformation, `k` 10, and `w` 0 (no flat bottom). The links of the last monomer of a polymer are ignored, unless it is `cyclic`. |
| `rotatable <a> <b>` | A bond from atom `a` of the monomer to atom `b` of the next one (indices local to each monomer, from 1, as in `link`), in place of links: it stands for the links that keep the length and the angles of the bond (the bond, and the bonds of `a` and `b` to their neighbors in their monomers), such that the monomers turn about it, as the [Rotatable bonds](@ref) of molecules. It requires the bonds of the polymer: the CONECT records of its file, or `bond_distance`. |
| `bond_distance <d>` | The bonds of the polymer, for `rotatable`: the pairs of atoms closer than `d` Å (in the first conformation; not pairs of hydrogens), if the file has no CONECT records. Without rotatable bonds, it has no effect. |
| `cyclic` | The last monomer is linked to the first one, with the links of the last monomer. |

The links of a junction keep its geometry, and leave free the rotation about the bond: as
many as needed (three, here, from the last two atoms of a monomer to the first two of
the next one, for coarse-grained beads). When the monomers are joined by a bond, as in
all-atom polymers, `rotatable` gives the links of the bond, from the bonds of the polymer:
five for a bond between carbon atoms with two other neighbors each. For example, for
polystyrene, with every bond of its backbone rotatable (see [Polymer-capped gold
nanoparticle](@ref)):

```
  bond_distance 1.9   # the bonds: the atoms closer than 1.9 Å
  repeat 20
    monomer 13        # CH(C6H5)
      rotatable 1 1   # the bond of atom 1 to atom 1 of the next monomer
    monomer 3         # CH2
      rotatable 1 1
  end repeat
```

- **The structure of the polymer.** The monomers have the geometry of each conformation
  of the file: the number of polymers is divided among the conformations, as for other
  molecules (see above), and each polymer is built from the monomers of its own.
- **Packing.** The polymers are first packed as rigid molecules, in their conformations,
  and then their monomers are released and the packing goes on with the flexible chains.
  (Chains released overlapping each other end up threaded through each other's monomers,
  and do not separate.) Polymers too dense to pack as rigid molecules are packed at a
  lower density, and compressed (see [Packing by compression](@ref)).
- **Conformations.** Any conformation of the polymer can be given (one is enough), as
  long as it fits, as a rigid molecule, in the region of the polymers: smaller than the
  periodic cell, and than the region of its constraints (expanded, with `use_compression`),
  since the polymers are first packed as rigid molecules. The flexibility of the chains
  does the rest. For example, a single extended conformation of the chains of the
  [Biomolecular condensate](@ref) (304 Å from end to end) packs in their homogeneous
  mixture, in a periodic cube of 260 Å (the chains fit along its diagonals), in a few
  seconds; but not in the droplet, a sphere of 100 Å of radius: rigid, the chains do
  not fit in it, and released into it from outside they fold into each other and jam.
  Nor does a larger expansion help: packed at 5% of the density (`use_compression 0.05`, in a
  sphere where the extended chains fit), they were compressed only to 41% of it. Coiled
  conformations, of the size of the chains in solution (e.g. from a simulation of the
  free chain), pack in both, the droplet compressed in less than a minute.
  Independently of the packing, several conformations are worth giving for the
  simulation: the polymers start from different conformations, and the results depend
  less on the choice of a single one (the flexibility of the packing only adapts the
  chains to their neighbors, and the simulation may take long to forget a conformation).
- **Close atoms.** The pairs of atoms of different monomers (up to three monomers apart)
  closer than the tolerance in a conformation (the bonded atoms of a junction, and their
  neighbors) are left out of the distances of the polymers of that conformation. With
  `rotatable` bonds, whose bonds are known, the pairs of atoms of different monomers up to
  three bonds apart are left out too.
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
  boundary conditions, a polymer is not wrapped across the boundaries of the cell: it is
  whole, placed in the cell by its center (the mean of the centers of its monomers), so
  large polymers cross the faces of the cell.
- **Not supported**: fixed polymers, `restart_from`, `constrain_rotation` and
  `voxel_place`.

Polymers are a case of flexible molecules, as molecules with [Rotatable bonds](@ref):
the monomers are the rigid fragments, given by their numbers of atoms, and the links are
written explicitly (the structure files of coarse-grained polymers rarely have bonds).

In Julia, the polymers are given by [`polymer_structure`](@ref), with the sequence made of
[`monomer`](@ref)s. The tutorial [Biomolecular condensate](@ref) packs a condensate of
polypeptides, in a homogeneous mixture and as a droplet, and the tutorial
[Polymer-capped gold nanoparticle](@ref) packs all-atom chains of polystyrene, flexible at
every bond of their backbone, bound to a fixed gold core.

## Rotatable bonds

Any molecule can be made flexible about some of its bonds: with
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
| `rotatable <i> <j> [<i> <j> ...]` | Rotatable bonds: pairs of bonded atoms, on one or several lines. The indices are the positions of the atoms in the structure file, from 1 (the first atom of the file is 1, whatever its serial number). A bond in a ring cannot rotate (an error). |
| `rotatable all` | Every rotatable bond of the molecule, found from its bonds (see below; for all-atom models with all their hydrogens). Other `rotatable` lines can add bonds to them. |
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

| Heads at (Å) | Rigid | 5 rotatable bonds | 28 rotatable bonds | 28, with `use_compression 0.5` |
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
and the heads are 23 to 29 Å from the ends of the chains. With `use_compression 0.5` (see
[Packing by compression](@ref)), the membrane is first packed 26% thicker, where the
rigid lipids fit, and then thinned with the lipids flexible: the rigid lipids, compressed
alone, stop at 82% of the density.

### All the rotatable bonds

With `rotatable all`, the rotatable bonds are found from the bonds of the molecule: the
bonds not in rings, between atoms that are not hydrogens and have each another bond to
an atom that is not a hydrogen (the bonds that turn more than hydrogen atoms), and that
are single for sure. The order of a bond is not in the structure file: it is deduced
from the numbers of bonds of its atoms. A carbon or nitrogen atom with fewer than four
bonds is, or may be, unsaturated (an amine nitrogen has three bonds, as an amide
nitrogen), and a bond between two such atoms does not turn: double and triple bonds,
amide bonds, and the bonds that are single but conjugated or might be double (aryl–vinyl,
aryl–carbonyl). The exception is the bonds between aromatic rings (rings of such atoms,
as in biphenyl), which turn. The bonds to saturated atoms turn: to carbon atoms with four
bonds, and to oxygen, sulfur or phosphorus atoms (ethers, esters, phosphates).

```
structure lignin.pdb
  number 7
  rotatable all
  inside sphere 0. 0. 0. 34.
end structure
```

!!! warning "All the hydrogen atoms are required"
    The rotatable bonds are found from the numbers of bonds of the atoms, which tell
    whether they are saturated only if the structure has all its hydrogen atoms. In a
    united-atom or coarse-grained model, or in a structure with some hydrogens missing,
    a CH₂ carbon has two bonds, and would be taken as unsaturated (an error if the
    molecule has no hydrogens; otherwise, check the bonds found). The number of bonds
    found is printed for each molecule. For other models, give the bonds by hand.

The bonds are those of the `CONECT` records of the file, or from `bond_distance`. For
POPC (with `bond_distance 1.8`), `rotatable all` finds 41 bonds: those of the chains,
the glycerol, the phosphate and the choline, but the double bond of the oleoyl chain. In
the oligomers of the [Lignin particle in acetone and water](@ref) tutorial, 24 to 72 bonds each.

In Julia, the molecules with rotatable bonds are given by [`flexible_structure`](@ref)
(`rotatable = :all` for `rotatable all`).
