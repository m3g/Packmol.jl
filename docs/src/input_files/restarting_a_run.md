# Restarting a run

`restart_from`/`restart_to` skip (or shortcut) the initial placement step and
resume packing from a previously-saved configuration — useful to extend a run
that stopped short of `nloop`, or to continue packing after adding more
molecules to the input file. As a global keyword, they cover every free
molecule in the system (and `restart_from` then skips the whole
initial-placement pipeline entirely, not just the random step — the main
reason to use it on a large system); given inside a `structure ... end
structure` block, they cover only that type's molecules instead, with the
rest of the system still placed normally.

Restarts use coordinate files (PDB, XYZ or mmCIF, by the extension of
`<file>`). `restart_to <file>` writes the coordinates of the system (or of the
structure type) each time the output is written, and `restart_from <file>`
reads any coordinate file with the same number of atoms and molecules as the
structure type(s) being restarted, in the same order: Packmol's own output, a
`restart_to` file, or a configuration from a simulation. Each molecule's
(center of mass, orientation) is recovered by rigid-body-aligning that
structure type's own template onto the observed atom positions (the Kabsch
algorithm), so the atoms don't need to sit at exactly the template's bond
lengths and angles.

The restart files of Fortran Packmol (centers of mass and Euler angles) are not
supported: restart from the coordinates of the previous run instead.

A whole-system restart file includes the fixed molecules, but their positions
are never read back: they are determined by their own `fixed` keyword.
