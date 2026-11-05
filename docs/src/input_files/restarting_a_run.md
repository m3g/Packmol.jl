# Restarting a run

`restart_from <file>` skips (or shortcuts) the initial placement step and
resumes packing from a previously saved configuration: the output of a previous
run, to extend a run that stopped short of `nloop` or to continue packing after
adding more molecules to the input file, or a configuration from a simulation. As
a global keyword, it covers every free molecule in the system (and then skips the
whole initial-placement pipeline, not just the random step: the main reason to use
it on a large system); given inside a `structure ... end structure` block, it covers
only that type's molecules, with the rest of the system still placed normally.

The file is a coordinate file (PDB, XYZ or mmCIF, by its extension), with the same
number of atoms and molecules as the structure type(s) being restarted, in the same
order. Each molecule's center of mass and orientation are recovered by
rigid-body-aligning the template of its structure type onto the observed atom
positions (the Kabsch algorithm), so the atoms don't need to sit at exactly the
template's bond lengths and angles. To restart a single structure type from the
output of a whole system, extract its atoms to a file first (with PDBTools, for
example).

The restart files of Fortran Packmol (centers of mass and Euler angles) are not
supported, and neither is its `restart_to` keyword (it is ignored, with a warning):
the output file itself is the restart file.

The output of a whole system includes the fixed molecules, but their positions
are never read back: they are determined by their own `fixed` keyword.
