# io

Input and output files.

| File | Contents |
|---|---|
| `read_input.jl` | The input file parser (`read_packmol_input`): the global keywords (and the deprecated and legacy ones of Fortran Packmol), and the structure blocks. |
| `read_structure_block.jl` | The contents of a structure block (`read_structure_data`): the structure file, the number of molecules, the constraints of the molecule and of `atoms` blocks, `constrain_rotation` blocks, and the other keywords. |
| `write_output.jl` | The packed atoms (`get_atoms`, with the molecules wrapped rigidly into the periodic cell) and the output file (`write_output`: PDB, XYZ or mmCIF, with the cell, `TER` and `CONECT` records, residue numbers and chains as Fortran Packmol). |
| `restart.jl` | Restarts (`restart_from`): the positions of the molecules recovered from coordinate files (PDB, XYZ or mmCIF) by aligning each molecule to its atoms. |

Structure files are read and written with PDBToolsBase.jl (`read_pdb`, `write_pdb`,
`read_mmcif`, ...), which is compatible with `juliac --trim`.
