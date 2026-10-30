#
# Structure files of the materials tutorials (docs/src/tutorials/files): graphite for
# the water drop, and the chabazite zeolite (from the IZA database, expanded to the full
# unit cell with Open Babel: `obabel CHA.cif -O CHA_uc.pdb --fillUC strict`) with methane
# and carbon dioxide.
#
using PDBTools
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files")

atom(name, resname, resnum, x, y, z; element=name[1:1]) =
    Atom(; name, resname, resnum, chain="A", x, y, z, pdb_element=element)

#
# Graphite: graphene layers in AB stacking, 3.35 Å apart, the top one at z = 0. The
# rectangular cell of graphene (2.46 x 4.26 Å) has 4 atoms.
#
function build_graphite(; nx=40, ny=23, nlayers=3, a=2.46, d=1.42)
    cell = [(0.0, 0.0), (a / 2, d / 2), (a / 2, 3d / 2), (0.0, 2d)]
    Lx, Ly = nx * a, ny * 3d
    atoms = Atom{Nothing}[]
    for layer in 0:nlayers-1, i in 0:nx-1, j in 0:ny-1, (x, y) in cell
        shift = isodd(layer) ? d : 0.0 # AB stacking
        push!(atoms, atom("C", "GRA", layer + 1,
            i * a + x - Lx / 2, mod(j * 3d + y + shift, Ly) - Ly / 2, -3.35 * layer))
    end
    return atoms
end

mkpath(joinpath(files, "water_drop"))
write_pdb(joinpath(files, "water_drop", "graphite.pdb"), build_graphite())
cp(joinpath(files, "membrane", "water.pdb"), joinpath(files, "water_drop", "water.pdb"); force=true)

#
# Chabazite (CHA), in its hexagonal cell (a = b = 13.675, c = 14.767 Å, γ = 120°):
# 2 x 2 x 2 unit cells, centered at the origin (as the `unitcell` of input files), with
# the CRYST1 record of the supercell. The unit cell has 36 Si and 72 O atoms.
#
mkpath(joinpath(files, "zeolite"))
uc = read_pdb(joinpath(@__DIR__, "CHA_uc.pdb"))
cell = Matrix{Float64}(UnitCell(13.675, 13.675, 14.767, 90, 90, 120))
A, B, C = cell[:, 1], cell[:, 2], cell[:, 3]
shift = -(A + B + C) # the center of the 2 x 2 x 2 supercell
cha = Atom{Nothing}[]
for i in 0:1, j in 0:1, k in 0:1, at in uc
    el = element(at) == "Si" ? "Si" : "O"
    x = position(at) + i * A + j * B + k * C + shift
    push!(cha, atom(el, "CHA", 1, x[1], x[2], x[3]; element=el))
end
s = Structure(cha)
s.unitcell = UnitCell(2 * 13.675, 2 * 13.675, 2 * 14.767, 90, 90, 120)
write_pdb(joinpath(files, "zeolite", "cha.pdb"), s)

# Methane (C-H 1.09 Å, tetrahedral) and carbon dioxide (C=O 1.16 Å)
h = 1.09 / sqrt(3)
write_pdb(joinpath(files, "zeolite", "methane.pdb"), [
    atom("C", "CH4", 1, 0.0, 0.0, 0.0),
    atom("H1", "CH4", 1, h, h, h), atom("H2", "CH4", 1, -h, -h, h),
    atom("H3", "CH4", 1, -h, h, -h), atom("H4", "CH4", 1, h, -h, -h),
])
write_pdb(joinpath(files, "zeolite", "co2.pdb"), [
    atom("C", "CO2", 1, 0.0, 0.0, 0.0),
    atom("O1", "CO2", 1, 0.0, 0.0, 1.16), atom("O2", "CO2", 1, 0.0, 0.0, -1.16),
])
