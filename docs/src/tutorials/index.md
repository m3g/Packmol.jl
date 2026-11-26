# Tutorials

The tutorials build complete systems, from a box with a mixture of two molecules to an
exosome with 13.6 million atoms. The first ones are the classic examples of
[Packmol](https://doi.org/10.1002/jcc.10216) (the bilayer and the vesicle redone
with an all-atom lipid, at realistic densities, the vesicle as nanoliposomes); the others
use more of the features of Packmol.jl: fixed molecules, periodic cells of different
shapes, curved surfaces, crystals. To have an AI agent build a system for you, see
[Using your AI agent](ai_agent.md), and the [Lipid nanoparticle](lnp.md), designed with an
agent from an illustration.

| Tutorial | System |
|---|---|
| [Water and urea](mixture.md) | A mixture of two molecules in a box |
| [Water–chloroform interface](interface.md) | Two solvents in adjacent boxes and a fixed molecule at the interface |
| [Lipid bilayer](bilayer.md) | A POPC bilayer in water, with periodic boundary conditions: lipids oriented by constraints on two of their atoms |
| [Nanoliposomes](nanoliposome.md) | Two POPC nanoliposomes in water, in a truncated octahedron periodic cell: lipid layers in spherical shells, and how curvature limits how many lipids fit in each leaflet |
| [Membrane tube with a neck](membrane_neck.md) | A POPC tube with a constriction, periodic along its axis: a shape given by a function, with every layer a custom constraint on the distance to the midplane of the bilayer |
| [Solvated protein](solvprotein.md) | A fixed protein in water and ions, in a rhombic dodecahedron periodic cell |
| [Protein in glycerol and water](protein_glycerol.md) | A fixed protein in a mixture, in a periodic box |
| [Water drop on graphite](water_drop.md) | A fixed surface, and a drop shaped by constraints |
| [Polymer in DMF](polymer_dmf.md) | A fixed polymer in a solvent, in a periodic box |
| [Carbon nanotube filled with water](nanotube.md) | A fixed nanotube crossing the periodic boundary, with water inside and outside |
| [Gases in a zeolite](zeolite.md) | A fixed crystal in its triclinic periodic cell, with gases in its pores |
| [ZIF-8 in water and ethanol](zif8.md) | A slab of a metal–organic framework in a liquid: ethanol in its pores, water kept out by a constraint; the numbers from `free_volume` |
| [Ionic liquid in a supercapacitor](supercapacitor.md) | An ionic liquid between charged graphite electrodes: the electric double layers set by the numbers of ions in regions |
| [Gold nanoparticle in toluene](gold_nanoparticle.md) | A fixed metal core with a shell of ligands, bound by one atom and oriented by another (constraints on atoms), in a solvent |
| [Polymer-capped gold nanoparticle](gold_polymer.md) | Flexible all-atom polymer chains (polystyrene) bound to a fixed gold core, in toluene, packed by compression around the fixed core |
| [Lignin particle in acetone and water](lignin.md) | A dense particle of flexible all-atom lignin oligomers (every rotatable bond found by `rotatable all`), packed by compression, then solvated around it; and dispersed in the solvent |
| [Liquid crystals](liquid_crystal.md) | Oriented molecules (rotation constraints): a nematic phase and a cholesteric (twisted) phase |
| [Biomolecular condensate](condensate.md) | Flexible chains (polymers) of a disordered protein in water, coarse-grained: a homogeneous mixture, and a droplet at 400 mg/mL packed by compression |
| [Zika virion](zika.md) | A large system: 2.2 million particles in membranes and spherical shells |
| [Lipid nanoparticle](lnp.md) | A system designed with an AI agent from an illustration: 5.6 million particles in monolayers, a lipid phase and water compartments with RNA |
| [Exosome](exosome.md) | A very large system, with all atoms: 13.6 million atoms in an asymmetric bilayer with membrane proteins, macromolecules inside, water and ions |
| [Minimal cell](cell.md) | The minimal cell JCVI-Syn3A of Bentopy (coarse-grained, about 40 million beads): 41,967 proteins, RNAs and metabolites placed by the voxel placement of Bentopy (`voxel_place`), then solvated |

## How the tutorials are organized

Each tutorial starts with a view of the final system and the files it needs. Then come
sections on the physics of the system: its structures, the regions of its molecules, and,
when the numbers of molecules are not trivial, how they follow from the densities and the
volumes of the regions. The sections that build the system (a drop, a phase, a step of a
large system) do it with an input file, run with the `packmol` executable, as with
Fortran Packmol (see [Installation](../installation.md)), and show how to check the
initial placement first.

The last section, **Julia script**, builds the same system with Julia code, with
[`structure_type`](@ref), [`PackmolSystem`](@ref) and [`packmol`](@ref) (see
[Defining systems in Julia](../julia_api/index.md)).

The input files and the Julia scripts are equivalent, and either can be followed alone. The input files of the
[Input files](../input_files/index.md) section and the functions of the
[Julia API](../julia_api/index.md) section are linked where they are first used.

## Getting the files

Each tutorial links to the files it needs: the input file and the structures of the
molecules (the Zika virion and the exosome, which are larger, are separate downloads). Put them in a
directory, and run the examples from it: the paths in an input file are relative to the
directory of the input file, and those in a Julia script to the current directory.

## [Checking the initial placement](@id check_initial_placement)

!!! tip "First check the initial placement of the molecules"
    Every input file of the tutorials has the line

    ```
    # check_initial_placement
    ```

    Uncomment it (remove the `#`), and run the input file: Packmol then only places each
    molecule in its region, satisfying the constraints, writes this initial placement to
    the output file, and stops, without packing. This takes a fraction of the time of the
    packing (2.5 minutes, for the exosome that takes 33 minutes). The molecules may still overlap,
    but a look at the output file in a molecular visualization program shows whether the
    constraints put them where intended. When they do, comment the line again, and run the
    packing.

    In Julia, the same is `check_initial_placement = true` in [`PackmolSystem`](@ref)
    (commented in the code of the tutorials).

    If the constraints cannot be satisfied at all, even without the distances between the
    molecules, Packmol stops after the initial placement with an error (unless the
    `continue_failed_initial_packing` keyword is set): that is almost always a sign of
    wrong constraints, cheaper to find before the packing than after it.

The tutorials of non-trivial systems show the initial placement before the packed system.
