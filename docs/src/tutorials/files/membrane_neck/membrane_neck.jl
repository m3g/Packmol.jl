#
# POPC membrane tube with a neck (a constriction), with water inside and outside, periodic
# along the axis of the tube (z). The shape is given by functions: the midplane of the
# bilayer is the surface of revolution r = ρ(z), and every layer is a custom constraint on
# the distance to it.
#
# Run from the directory of popc.pdb and water.pdb, with several threads:
#
#     julia -t 14 membrane_neck.jl
#
using Packmol

# Geometry (Å)
const R = 80.0     # radius of the midplane of the tube, far from the neck
const Rn = 48.0    # radius of the midplane at the neck
const σ = 35.0     # width of the neck
const L = 300.0    # length of the periodic cell, along the tube
const a = 2 * (R + 30 + 15) # side of the cell, perpendicular to the tube

# The midplane, r = ρ(z), and its slope
ρ(z) = R - (R - Rn) * exp(-z^2 / (2σ^2))
dρ(z) = (R - Rn) * z / σ^2 * exp(-z^2 / (2σ^2))

# Distance of a point to the midplane, positive outside the tube (to first order: the
# radial distance, corrected by the slope of the surface)
h(x) = (hypot(x[1], x[2]) - ρ(x[3])) / sqrt(1 + dρ(x[3])^2)

# Regions beyond (h ≥ d) and within (h ≤ d) a distance d of the midplane
beyond(d) = CustomConstraint(x -> d - h(x))
within(d) = CustomConstraint(x -> h(x) - d)

# Areas of the surfaces at a distance d of the midplane, and volumes inside them
const zs = range(-L/2, L/2; length = 3001)
area(d) = sum(2π * (ρ(z) + d) * sqrt(1 + dρ(z)^2) for z in zs) * step(zs)
volume(d) = sum(π * (ρ(z) + d)^2 for z in zs) * step(zs)

# Lipids: 64 Å² per lipid at the narrow surface of each leaflet (the inner heads, and the
# midplane for the outer leaflet), plus 10%. Water at 1 g/mL: 29.9 Å³ per molecule.
n_inner = round(Int, 1.1 * area(-28) / 64)
n_outer = round(Int, 1.1 * area(0) / 64)
n_water_inside = round(Int, volume(-31) / 29.9)
n_water_outside = round(Int, (a^2 * L - volume(30)) / 29.9)

molecules = [
    structure_type("water.pdb"; number = n_water_inside, constraints = [within(-31.0)]),
    # Inner leaflet: heads (atom 1) towards the axis, chain ends (88, 131) near the midplane
    structure_type("popc.pdb"; number = n_inner,
        atom_constraints = [[1] => within(-28.0), [88, 131] => beyond(-6.0)]),
    # Outer leaflet: chain ends near the midplane, heads outwards
    structure_type("popc.pdb"; number = n_outer,
        atom_constraints = [[88, 131] => within(6.0), [1] => beyond(28.0)]),
    structure_type("water.pdb"; number = n_water_outside, constraints = [beyond(30.0)]),
]
system = PackmolSystem(molecules;
    output = "membrane_neck.pdb", tolerance = 2.0,
    unitcell = [a, a, L, 90, 90, 90], # orthorhombic, centered at the origin
    # uncomment to only place the molecules
    # check_initial_placement = true,
)
packmol(system)
