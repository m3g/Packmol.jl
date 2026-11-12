# tools

Tools for setting up systems, used before packing.

| File | Contents |
|---|---|
| `volume.jl` | `guess_volume` (the volume of molecules at a density) and `free_volume` (the volume of the region of each structure type, by Monte Carlo, and the density its molecules would have); also the `--guess-volume` and `--free-volume` options of the command line. |
| `shape.jl` | `shape_descriptors` (asphericity, acylindricity and anisotropy, from the gyration tensor) and `select_conformations` (the most rod-like or sphere-like conformations of a multi-model file). |
| `voxel_place.jl` | `voxel_place`: a direct port to Julia of the placement algorithm of [Bentopy](https://github.com/marrink-lab/bentopy) (its packing step): molecules voxelized at random rotations, placed at random free voxels of their compartments (the constraints as masks), the largest moment first, those that do not fit dropped; returns the system of the placed molecules, which `packmol` starts from. |
| `system_building.jl` | Public, not exported: `builtin_lipid` (the files, head and chain-end atoms, area and length of the built-in lipids) and `cross_section` (the area of a membrane protein in a leaflet). |
