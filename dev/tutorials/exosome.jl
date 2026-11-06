#
# Builds the files of the exosome tutorial (the packmol_exosome set of MolSimTestFiles):
# an extracellular vesicle with tetraspanins (CD9 and CD81) in an asymmetric bilayer,
# Argonaute-2/miRNA complexes and nucleosomes inside, water and NaCl inside and outside,
# in a truncated octahedron. Two input files: `proteins.inp`, which packs the proteins
# (the tetraspanins oriented across the membrane, the complexes and the nucleosomes in
# the lumen), and `exosome.inp`, which packs the lipids, the water and the ions around
# them, fixed.
#
#     julia --project=dev/large_tests dev/tutorials/exosome.jl <prepared> <outdir>
#
# `<prepared>` has the protein and nucleic acid structures with hydrogens (mmCIF):
# cd81.cif and cd9.cif (OPM 5TCX and 6K4J: membrane normal along z, middle of the
# membrane at z = 0, extracellular side at z > 0), ago2_mirna.cif (PDB 4W5N) and
# nucleosome.cif (PDB 1KX5), prepared with PDBFixer (missing heavy atoms and hydrogens
# at pH 7, no missing residues built, waters, ions, ligands and DUM atoms removed).
#
using Packmol, StaticArrays, LinearAlgebra, Printf
import PDBToolsBase
using PDBToolsBase: read_pdb, read_mmcif, write_mmcif, mass

const prepared, outdir = ARGS[1], ARGS[2]
mkpath(outdir)

# Geometry (Å): the middle of the bilayer, and the periodic cell
const Rm = 228.0            # radius of the middle of the bilayer
const d = 580.0             # distance between opposite (hexagonal) faces of the cell
const chain_end_gap = 6.0   # chain ends within 6 Å of the middle of the bilayer
const flexibility = 0.15    # heads beyond chain_end_gap + (1 - flexibility) L
const solvent_gap = 2.0     # water beyond the heads (outside) / 3 Å (inside)
const tolerance = 2.0

# Lipids: outer leaflet rich in PC (DPPC standing for the saturated sphingomyelin), inner
# leaflet rich in PE and PS, as the plasma membrane from which exosomes bud
# (the built-in lipids: structure file, head and chain-end atoms, area, length)
lipids = Dict(name => Packmol.builtin_lipid(name) for name in ("POPC", "POPE", "POPS", "DPPC"))
outer = ["POPC" => 2, "DPPC" => 1]
inner = ["POPE" => 2, "POPS" => 1, "POPC" => 1]
lipid_charge = Dict("POPC" => 0, "POPE" => 0, "POPS" => -1, "DPPC" => 0)

const lipid_dir = joinpath(Packmol.RecipesDirectory, "lipids")
coords(atoms) = [SVector{3,Float64}(a.x, a.y, a.z) for a in atoms]
head_offset(name) = chain_end_gap + (1 - flexibility) * lipids[name].length

# Radii of the heads and of the water
r_inner_heads = Rm - maximum(head_offset(first(p)) for p in inner)
r_outer_heads = Rm + maximum(head_offset(first(p)) for p in outer)
r_water_in = r_inner_heads - 3.0
r_water_out = r_outer_heads + solvent_gap

# Membrane proteins: 40 copies of CD81 and of CD9, oriented across the membrane as the
# lipids: the Cα atoms 18 to 24 Å from the middle of the membrane on the extracellular
# side (z > 0 in OPM) outside a sphere of radius Rm + 17, and those on the intracellular
# side inside a sphere of radius Rm - 17 (the middle of their membrane within a few Å of Rm,
# tilted by 10° on average)
const n_cd81, n_cd9 = 40, 40
function membrane_protein(file)
    atoms = read_mmcif(joinpath(prepared, file))
    # PDBFixer adds an OXT to the last residue on top of its O (0.6 Å apart): removed
    last_residue = atoms[end].resnum
    atoms = [a for a in atoms if !(a.name == "OXT" && a.resnum == last_residue)]
    write_mmcif(joinpath(outdir, file), atoms)
    band(lo, hi) = [i for (i, a) in enumerate(atoms) if a.name == "CA" && lo <= a.z <= hi]
    return atoms, coords(atoms), band(18.0, 24.0), band(-24.0, -18.0)
end
cd81, x_cd81, out_cd81, in_cd81 = membrane_protein("cd81.cif")
cd9, x_cd9, out_cd9, in_cd9 = membrane_protein("cd9.cif")

# Lipids per leaflet: the safe count of the nanoliposome tutorial (the area of the most
# crowded surface of the leaflet: the heads, inside, and the middle of the bilayer,
# outside), plus 10%, without the cross-sections of the proteins in the leaflet

section(x, z1, z2) = Packmol.cross_section(x; z = (z1, z2), radius = tolerance)
# (in the hydrophobic core, 17 Å on each side of the middle: the half-thickness of OPM)
const core = 17.0
sections_out = n_cd81 * section(x_cd81, 0.0, core) + n_cd9 * section(x_cd9, 0.0, core)
sections_in = n_cd81 * section(x_cd81, -core, 0.0) + n_cd9 * section(x_cd9, -core, 0.0)
mean_area(composition) = sum(lipids[l].area * w for (l, w) in composition) / sum(last, composition)
n_outer = round(Int, 1.1 * (4π * Rm^2 - sections_out) / mean_area(outer))
n_inner = round(Int, 1.1 * (4π * r_inner_heads^2 - sections_in) / mean_area(inner))
split_counts(n, composition) = begin
    w = [last(p) for p in composition]
    c = round.(Int, n .* w ./ sum(w))
    c[1] += n - sum(c)
    [first(p) => ci for (p, ci) in zip(composition, c)]
end
outer_counts = split_counts(n_outer, outer)
inner_counts = split_counts(n_inner, inner)

# Contents of the lumen: in a sphere 5 Å inside the water of the lumen
const n_ago2, n_nucleosome = 8, 2
ago2 = read_mmcif(joinpath(prepared, "ago2_mirna.cif"))
nucleosome = read_mmcif(joinpath(prepared, "nucleosome.cif"))
write_mmcif(joinpath(outdir, "ago2_mirna.cif"), ago2)
write_mmcif(joinpath(outdir, "nucleosome.cif"), nucleosome)
r_lumen = r_water_in - 5

# Charges of the residues at pH 7 (termini not counted; nucleotides -1)
function charge(atoms)
    q = 0
    for r in unique((a.chain, a.resnum, a.resname) for a in atoms)
        q += r[3] in ("ARG", "LYS") ? 1 : r[3] in ("ASP", "GLU") ? -1 :
             r[3] in ("A", "U", "G", "C", "DA", "DT", "DG", "DC") ? -1 : 0
    end
    return q
end
q_cd81, q_cd9, q_ago2, q_nucleosome = charge(cd81), charge(cd9), charge(ago2), charge(nucleosome)

# Volumes (Å³): of the macromolecules, at 1.35 g/mL; of the water inside and outside
volume(atoms) = mass(atoms) / (0.60221 * 1.35)
cell_volume = abs(det(octahedral_unitcell(Float64, d)))
outside_fraction(x) = count(p -> p[3] > r_outer_heads - Rm, x) / length(x)
v_proteins_out = n_cd81 * volume(cd81) * outside_fraction(x_cd81) + n_cd9 * volume(cd9) * outside_fraction(x_cd9)
v_in = 4π / 3 * r_water_in^3 - n_ago2 * volume(ago2) - n_nucleosome * volume(nucleosome)
v_out = cell_volume - 4π / 3 * r_water_out^3 - v_proteins_out

# Water at 1 g/mL (33.4 molecules/nm³), NaCl at 0.15 mol/L, and counterions: inside, for
# the contents of the lumen and the PS of the inner leaflet; outside, for the proteins
water_density = 0.0334    # molecules / Å³
salt = 0.15 * 6.02214e-4  # ion pairs / Å³
n_pops = sum(c for (l, c) in inner_counts if l == "POPS"; init = 0)
q_in = n_ago2 * q_ago2 + n_nucleosome * q_nucleosome + n_pops * lipid_charge["POPS"]
q_out = n_cd81 * q_cd81 + n_cd9 * q_cd9
ions(v, q) = (round(Int, salt * v) + max(0, -q), round(Int, salt * v) + max(0, q))
na_in, cl_in = ions(v_in, q_in)
na_out, cl_out = ions(v_out, q_out)
water_in = round(Int, water_density * v_in) - na_in - cl_in
water_out = round(Int, water_density * v_out) - na_out - cl_out

# Small structures, copied as mmCIF
for name in ("water", "sodium", "chloride", "popc", "pope", "pops", "dppc")
    write_mmcif(joinpath(outdir, "$name.cif"), read_pdb(joinpath(lipid_dir, "$name.pdb")))
end

# The input file
f2(x) = @sprintf("%.1f", x)

# Step 1: the proteins
membrane_block(file, n, out, in) = """
structure $file
  number $n
  atoms $(join(out, " "))
    outside sphere 0. 0. 0. $(f2(Rm + 17))
  end atoms
  atoms $(join(in, " "))
    inside sphere 0. 0. 0. $(f2(Rm - 17))
  end atoms
end structure
"""
write(joinpath(outdir, "proteins.inp"), """
#
# The proteins of the exosome: 40 copies of each of the tetraspanins CD81 and CD9,
# crossing the membrane (the middle of the membrane at $(f2(Rm)) Å from the center), and,
# in the lumen, 8 complexes of Argonaute-2 with a miRNA and 2 nucleosomes. The output
# is fixed in the packing of the exosome (exosome.inp).
#
tolerance $tolerance
output proteins.cif

# To write only the initial placement of the molecules to the output file,
# without packing them (to check the constraints), uncomment the line below.

# check_initial_placement

# The tetraspanins: the Cα atoms 18 to 24 Å from the middle of the membrane on the
# extracellular side (first atoms block) outside, and on the intracellular side (second
# atoms block) inside
$(membrane_block("cd81.cif", n_cd81, out_cd81, in_cd81))
$(membrane_block("cd9.cif", n_cd9, out_cd9, in_cd9))
# The lumen
structure ago2_mirna.cif
  number $n_ago2
  inside sphere 0. 0. 0. $(f2(r_lumen))
end structure

structure nucleosome.cif
  number $n_nucleosome
  inside sphere 0. 0. 0. $(f2(r_lumen))
end structure
""")

# Step 2: the lipids, water and ions around the proteins
io = IOBuffer()
println(io, """
#
# An exosome: an extracellular vesicle of 51 nm, with 40 copies of each of the
# tetraspanins CD81 and CD9 in an asymmetric bilayer, 8 complexes of Argonaute-2 with
# a miRNA and 2 nucleosomes inside, and water and NaCl (0.15 mol/L) inside and outside,
# in a truncated octahedron. All atoms, written as mmCIF.
#
tolerance $tolerance
output exosome.cif
pbc octahedral 0. 0. 0. $(f2(d))
periodic_boundary_style octahedral

# To write only the initial placement of the molecules to the output file,
# without packing them (to check the constraints), uncomment the line below.

# check_initial_placement

# The proteins, packed by proteins.inp, fixed
structure proteins.cif
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure
""")
for (l, c) in inner_counts
    p = lipids[l]
    println(io, """
    # Inner leaflet: $l
    structure $(lowercase(l)).cif
      number $c
      atoms $(p.head)
        inside sphere 0. 0. 0. $(f2(Rm - head_offset(l)))
      end atoms
      atoms $(join(p.chain_ends, " "))
        outside sphere 0. 0. 0. $(f2(Rm - chain_end_gap))
      end atoms
    end structure
    """)
end
for (l, c) in outer_counts
    p = lipids[l]
    println(io, """
    # Outer leaflet: $l
    structure $(lowercase(l)).cif
      number $c
      atoms $(join(p.chain_ends, " "))
        inside sphere 0. 0. 0. $(f2(Rm + chain_end_gap))
      end atoms
      atoms $(p.head)
        outside sphere 0. 0. 0. $(f2(Rm + head_offset(l)))
      end atoms
    end structure
    """)
end
println(io, """
# The lumen: water and ions
structure water.cif
  number $water_in
  inside sphere 0. 0. 0. $(f2(r_water_in))
end structure

structure sodium.cif
  number $na_in
  inside sphere 0. 0. 0. $(f2(r_water_in))
end structure

structure chloride.cif
  number $cl_in
  inside sphere 0. 0. 0. $(f2(r_water_in))
end structure

# Outside: water and ions
structure water.cif
  number $water_out
  outside sphere 0. 0. 0. $(f2(r_water_out))
end structure

structure sodium.cif
  number $na_out
  outside sphere 0. 0. 0. $(f2(r_water_out))
end structure

structure chloride.cif
  number $cl_out
  outside sphere 0. 0. 0. $(f2(r_water_out))
end structure""")
write(joinpath(outdir, "exosome.inp"), String(take!(io)))

natoms = n_cd81 * length(cd81) + n_cd9 * length(cd9) + n_ago2 * length(ago2) + n_nucleosome * length(nucleosome) +
    3 * (water_in + water_out) + na_in + na_out + cl_in + cl_out +
    sum(c * length(read_pdb(joinpath(lipid_dir, lowercase(l) * ".pdb"))) for (l, c) in [inner_counts; outer_counts])
println("lipid lengths: ", Dict(k => round(v.length, digits = 1) for (k, v) in lipids))
println("radii: inner heads < ", f2(r_inner_heads), ", outer heads > ", f2(r_outer_heads), ", water in < ", f2(r_water_in), ", out > ", f2(r_water_out))
println("protein sections: outer ", round(Int, sections_out), " Å², inner ", round(Int, sections_in), " Å²")
println("lipids: inner ", inner_counts, ", outer ", outer_counts)
println("charges: CD81 $q_cd81, CD9 $q_cd9, Ago2-miRNA $q_ago2, nucleosome $q_nucleosome; inside $q_in, outside $q_out")
println("water: inside $water_in, outside $water_out; Na+ $na_in/$na_out, Cl- $cl_in/$cl_out")
println("total: ", natoms, " atoms; cell volume ", round(cell_volume / 1e6, digits = 1), " nm³ × 10³")
