#
# MembraneProtein: a membrane protein in a lipid bilayer, solvated by water and ions,
# under orthorhombic periodic boundary conditions. The protein is given oriented for the
# membrane, as in the OPM database: the membrane normal along z, and the middle of the
# bilayer at z = 0. It is kept fixed, and the lipids, placed as in `Membrane`, fill the
# rest of the membrane; the water and the ions fill the slabs beyond the heads of the
# lipids.
#
export MembraneProtein

mutable struct MembraneProtein <: Recipe
    protein_pdbfile::String
    protein_charge::Int
    membrane::Membrane # the lipids (composition, area per lipid, flexibility)
    margin::typeof(1.0u"Å")
    solvent_margin::typeof(1.0u"Å")
    ionic_concentration::typeof(1.0u"mol/L")
    water_pdbfile::String
    cation_pdbfile::String
    anion_pdbfile::String
    cation_charge::Int
    anion_charge::Int
    water_molar_mass::typeof(1.0u"g/mol")
    cation_molar_mass::typeof(1.0u"g/mol")
    anion_molar_mass::typeof(1.0u"g/mol")
    density_table::DensityTable
end

_mandatory_keywords(::Type{MembraneProtein}) = (:protein, :lipids)

"""
    MembraneProtein(;
        protein::String,
        lipids::Union{String,Vector{String}},
        lipid_molar_ratio = nothing, # equal amounts if not given
        margin::Number = 15.0, # Å, between the protein and the faces of the cell, in x and y
        solvent_margin::Number = 15.0, # Å, of water beyond the protein and the heads, in z
        ionic_concentration::Number = 0.15, # mol/L
        protein_charge = nothing, # computed from the residues of the protein by default
        area_per_lipid = nothing, # Å², automatic for the built-in lipids
        lipid_head = nothing, lipid_tail = nothing, lipid_molar_mass = nothing,
        flexibility::Real = 0.15,
        water_pdbfile = nothing, # water by default
        cation_pdbfile = nothing, anion_pdbfile = nothing, # Na⁺ and Cl⁻ by default
        cation_charge = nothing, anion_charge = nothing,
        density_table = nothing, density_table_units = "mol/kg", # NaCl(aq) by default
    )

Setup a membrane protein in a lipid bilayer, solvated by water and ions, with orthorhombic
periodic boundary conditions. The `protein` (a PDB file) must be oriented for the membrane,
as the structures of the [OPM database](https://opm.phar.umich.edu) are: with the membrane
normal along z, and the middle of the bilayer at z = 0 (remove the `DUM` atoms of the OPM
files, which mark the boundaries of the membrane). The protein is kept fixed, at the
coordinates of the file.

The lipids are those of [`Membrane`](@ref): the names of built-in lipids (`"POPC"`,
`"POPE"`, `"POPS"`, `"POPG"`, `"DPPC"`, `"DOPC"`, `"DMPC"`) or PDB files, mixed at
`lipid_molar_ratio`, and placed in the same way. The cell extends `margin` beyond the
protein in x and y, and the number of lipids in each leaflet fills its area, without the
cross-section of the protein in the leaflet, at `area_per_lipid`.

The water and the ions fill two slabs, beyond the heads of the lipids of each leaflet; the
cell extends `solvent_margin` beyond the heads and beyond the protein, in z. The water and
the ions are placed in the volume of the slabs without the protein (estimated by Monte
Carlo, as [`free_volume`](@ref)): ions at `ionic_concentration` (mol/L), and counterions to
neutralize the charge of the protein, and water to the density of the salt solution, as in
[`SolutionBoxUWI`](@ref).
"""
function MembraneProtein(;
    protein::String = _missing_keyword(MembraneProtein, :protein),
    lipids::Union{String,Vector{String}} = _missing_keyword(MembraneProtein, :lipids),
    lipid_molar_ratio::Union{Nothing,Vector{<:Real}} = nothing,
    margin::Number = 15.0,
    solvent_margin::Number = 15.0,
    ionic_concentration::Number = 0.15,
    protein_charge::Union{Nothing,Integer} = nothing,
    area_per_lipid::Union{Nothing,Number} = nothing,
    lipid_head::Union{Nothing,Vector{<:Integer}} = nothing,
    lipid_tail::Union{Nothing,Vector} = nothing,
    lipid_molar_mass::Union{Nothing,Vector{<:Number}} = nothing,
    flexibility::Real = 0.15,
    water_pdbfile::Union{Nothing,String} = nothing,
    cation_pdbfile::Union{Nothing,String} = nothing,
    anion_pdbfile::Union{Nothing,String} = nothing,
    cation_charge::Union{Nothing,Integer} = nothing,
    anion_charge::Union{Nothing,Integer} = nothing,
    density_table::Union{Nothing,Matrix{<:Number}} = nothing,
    density_table_units::String = "mol/kg",
)
    protein = expanduser(protein)
    isfile(protein) || throw(ArgumentError("Protein file $protein not found."))
    protein_atoms = read_pdb(protein)
    isnothing(protein_charge) && (protein_charge = round(Int, sum(charge, eachresidue(protein_atoms))))

    # The lipids (the size of this membrane is not used)
    membrane = Membrane(; lipids, lipid_molar_ratio, lipid_head, lipid_tail, area_per_lipid,
        lipid_molar_mass, flexibility, total_lipids = 2)

    margin = _ensure_unit(margin, u"Å")
    solvent_margin = _ensure_unit(solvent_margin, u"Å")
    ustrip(margin) >= 0 || throw(ArgumentError("margin must be non-negative."))
    ustrip(solvent_margin) > 0 || throw(ArgumentError("solvent_margin must be positive."))
    unit(ionic_concentration) == NoUnits && (ionic_concentration = ionic_concentration * u"mol/L")
    ionic_concentration = _ensure_unit(ionic_concentration, u"mol/L")
    ustrip(ionic_concentration) >= 0 || throw(ArgumentError("ionic_concentration must be non-negative."))

    # Water and ions: built-in by default (the charges of custom ions must be given)
    custom_ions = !isnothing(cation_pdbfile) || !isnothing(anion_pdbfile)
    water_pdbfile = isnothing(water_pdbfile) ? joinpath(_MEMBRANE_LIPIDS_DIR, "water.pdb") : expanduser(water_pdbfile)
    if isnothing(cation_pdbfile)
        cation_pdbfile = joinpath(_MEMBRANE_LIPIDS_DIR, "sodium.pdb")
        isnothing(cation_charge) && (cation_charge = 1)
    else
        cation_pdbfile = expanduser(cation_pdbfile)
        isnothing(cation_charge) && throw(ArgumentError("cation_charge must be given with cation_pdbfile."))
    end
    if isnothing(anion_pdbfile)
        anion_pdbfile = joinpath(_MEMBRANE_LIPIDS_DIR, "chloride.pdb")
        isnothing(anion_charge) && (anion_charge = -1)
    else
        anion_pdbfile = expanduser(anion_pdbfile)
        isnothing(anion_charge) && throw(ArgumentError("anion_charge must be given with anion_pdbfile."))
    end
    cation_charge > 0 || throw(ArgumentError("cation_charge must be positive, got $cation_charge."))
    anion_charge < 0 || throw(ArgumentError("anion_charge must be negative, got $anion_charge."))
    water_molar_mass = _ensure_unit(mass(read_pdb(water_pdbfile)), u"g/mol")
    cation_molar_mass = _ensure_unit(mass(read_pdb(cation_pdbfile)), u"g/mol")
    anion_molar_mass = _ensure_unit(mass(read_pdb(anion_pdbfile)), u"g/mol")
    if isnothing(density_table)
        custom_ions && @warn "Custom ions without a `density_table`: the density of aqueous NaCl is assumed." _file=nothing _line=nothing
        density_table = density_NaCl_aq
    end
    dtable = DensityTable(density_table_units, density_table[:, 1], density_table[:, 2],
        water_molar_mass, cation_molar_mass + anion_molar_mass)

    return MembraneProtein(protein, protein_charge, membrane, margin, solvent_margin, ionic_concentration,
        water_pdbfile, cation_pdbfile, anion_pdbfile, cation_charge, anion_charge,
        water_molar_mass, cation_molar_mass, anion_molar_mass, dtable)
end

function Base.show(io::IO, ::MIME"text/plain", system::MembraneProtein)
    m = system.membrane
    print(io, chomp("""
    ==================================================================
    Membrane protein:
    ==================================================================
        Protein: $(basename(system.protein_pdbfile)) (charge $(system.protein_charge))
        Lipids: $(join(m.lipid_names, ", "))
        Lipid molar ratio (normalized): $(round.(m.lipid_weight; digits=4))
        Area per lipid: $(round(u"Å^2", m.area_per_lipid; digits=2))
        Margin around the protein: $(system.margin)
        Solvent margin: $(system.solvent_margin)
        Ionic concentration: $(system.ionic_concentration)
        Ions: $(basename(system.cation_pdbfile)) (+$(system.cation_charge)), $(basename(system.anion_pdbfile)) ($(system.anion_charge))
    ==================================================================
    """))
end

# The area (Å²) of the cross-section of the atoms `xs` in the rectangle [lo, hi] of the
# xy plane: the union of discs of radius `r` around them, on a grid of spacing `h`
function _cross_section(xs::Vector{SVector{3,Float64}}, lo, hi, r; h=0.25)
    nx, ny = ceil(Int, (hi[1] - lo[1]) / h), ceil(Int, (hi[2] - lo[2]) / h)
    grid = falses(nx, ny)
    k = ceil(Int, r / h)
    for x in xs
        i0, j0 = floor(Int, (x[1] - lo[1]) / h) + 1, floor(Int, (x[2] - lo[2]) / h) + 1
        for i in i0-k:i0+k, j in j0-k:j0+k
            (1 <= i <= nx && 1 <= j <= ny) || continue
            cx, cy = lo[1] + (i - 0.5) * h, lo[2] + (j - 0.5) * h
            (cx - x[1])^2 + (cy - x[2])^2 <= r^2 && (grid[i, j] = true)
        end
    end
    return count(grid) * h^2
end

function _setup(system::MembraneProtein; npoints::Int = 200_000)
    (; protein_pdbfile, protein_charge, membrane, margin, solvent_margin, ionic_concentration,
       water_pdbfile, cation_pdbfile, anion_pdbfile, cation_charge, anion_charge,
       water_molar_mass, cation_molar_mass, anion_molar_mass, density_table) = system
    (; lipid_names, lipid_pdbfiles, lipid_head, lipid_tail, lipid_head_tail_length, lipid_weight,
       area_per_lipid, flexibility) = membrane
    tolerance = 2.0
    protein = read_pdb(protein_pdbfile)
    xs = [SVector{3,Float64}(a.x, a.y, a.z) for a in protein]
    plo = reduce((a, b) -> min.(a, b), xs)
    phi = reduce((a, b) -> max.(a, b), xs)

    # The planes of the lipids, from the middle of the bilayer (z = 0), as in `Membrane`
    nlip = length(lipid_pdbfiles)
    L = ustrip.(u"Å", lipid_head_tail_length)
    head_planes = [_CHAIN_END_PLANE + (1 - flexibility) * L[i] for i in 1:nlip]
    zsolvent = maximum(head_planes) + _SOLVENT_GAP

    # The cell: `margin` beyond the protein in x and y, `solvent_margin` beyond the heads
    # and the protein in z
    m, sm = ustrip(u"Å", margin), ustrip(u"Å", solvent_margin)
    lo = SVector(plo[1] - m, plo[2] - m, min(-zsolvent, plo[3]) - sm)
    hi = SVector(phi[1] + m, phi[2] + m, max(zsolvent, phi[3]) + sm)
    sides = hi - lo
    unitcell = Matrix{Float64}(Diagonal(collect(sides)))
    center = (lo + hi) / 2

    # Lipids per leaflet: the area of the cell without the cross-section of the protein in
    # the leaflet (its atoms within the tolerance), at the area per lipid
    area = sides[1] * sides[2]
    apl = ustrip(u"Å^2", area_per_lipid)
    leaflets = ((true, 0.0, zsolvent), (false, -zsolvent, 0.0))
    lipid_placements = _LipidPlacement[]
    leaflet_lines = String[]
    for (upper, z1, z2) in leaflets
        section = _cross_section([x for x in xs if z1 <= x[3] <= z2], lo, hi, tolerance)
        n = round(Int, (area - section) / apl)
        n > 0 || throw(ArgumentError("No lipids fit in the $(upper ? "upper" : "lower") leaflet: increase the margin."))
        counts = round.(Int, lipid_weight .* n)
        counts[end] += n - sum(counts)
        s = upper ? 1.0 : -1.0
        for i in 1:nlip
            counts[i] == 0 && continue
            push!(lipid_placements, _LipidPlacement(lipid_pdbfiles[i], counts[i], lipid_head[i], lipid_tail[i],
                s * head_planes[i], s * _CHAIN_END_PLANE, upper, nothing))
        end
        push!(leaflet_lines, "    $(upper ? "upper" : "lower") leaflet: $n lipids ($(join(["$(counts[i]) $(lipid_names[i])" for i in 1:nlip], ", "))), " *
            "protein cross-section $(round(section; digits=1)) Å²")
    end

    # The volume of the solvent slabs without the protein, by Monte Carlo
    slabs = ((zsolvent, hi[3]), (lo[3], -zsolvent))
    rng = Random.Xoshiro(1)
    slab_volumes = Float64[]
    uc = SMatrix{3,3,Float64}(unitcell)
    for (z1, z2) in slabs
        points = [SVector(lo[1] + rand(rng) * sides[1], lo[2] + rand(rng) * sides[2], z1 + rand(rng) * (z2 - z1)) for _ in 1:npoints]
        near = _near_fixed(points, xs, tolerance, true, uc)
        push!(slab_volumes, (npoints - count(near)) / npoints * area * (z2 - z1))
    end
    vs = sum(slab_volumes) * u"Å^3"

    # Ions and water, as in `SolutionBoxUWI`
    Msalt = cation_molar_mass + anion_molar_mass
    nsalt = round(Int, Unitful.Na * ionic_concentration * vs)
    ncation = nsalt * abs(anion_charge)
    nanion = nsalt * abs(cation_charge)
    if protein_charge > 0
        mod(protein_charge, abs(anion_charge)) == 0 || throw(ArgumentError(
            "The charge of the protein ($protein_charge) is not a multiple of that of the anion ($anion_charge)."))
        nanion += protein_charge ÷ abs(anion_charge)
    elseif protein_charge < 0
        mod(-protein_charge, cation_charge) == 0 || throw(ArgumentError(
            "The charge of the protein ($protein_charge) is not a multiple of that of the cation ($cation_charge)."))
        ncation += (-protein_charge) ÷ cation_charge
    end
    ρ = ionic_concentration <= first(density_table.concentration) ? PURE_WATER_DENSITY :
        interpolate_density(density_table, water_molar_mass, Msalt, ionic_concentration, "mol/L")
    mions = (ncation / Unitful.Na) * cation_molar_mass + (nanion / Unitful.Na) * anion_molar_mass
    mwater = uconvert(u"g", ρ * vs) - uconvert(u"g", mions)
    mwater > zero(mwater) || throw(ArgumentError("The mass of the ions exceeds that of the solution."))
    nwater = round(Int, Unitful.Na * mwater / water_molar_mass)

    # Water and ions in the two slabs, in proportion to their volumes
    f = slab_volumes[1] / sum(slab_volumes)
    split(n) = (k = round(Int, f * n); (k, n - k))
    solvent_placements = _SolventPlacement[]
    for (file, n) in ((water_pdbfile, nwater), (cation_pdbfile, ncation), (anion_pdbfile, nanion))
        for (k, (z1, z2)) in zip(split(n), slabs)
            k > 0 && push!(solvent_placements, _SolventPlacement(file, k, z1, z2))
        end
    end

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Protein: $(basename(protein_pdbfile)), charge $protein_charge
        Periodic box = $(_unitcell_description(unitcell))
        Center of the box = $(round.(center; digits=3))

        Lipids (area per lipid = $(round(apl; digits=2)) Å²):
        $(join(leaflet_lines, "\n"))

        Solvent ($(round(ustrip(u"Å^3", vs); digits=0)) Å³ without the protein):
            water: $nwater
            $(basename(cation_pdbfile)): $ncation
            $(basename(anion_pdbfile)): $nanion
            density of the solution: $(round(u"g/mL", ρ; digits=4))
            total charge: $(protein_charge + ncation * cation_charge + nanion * anion_charge)

        ==================================================================
        """
    return (; lipid_placements, solvent_placements, unitcell, center, lo, hi, summary)
end

"""
    write_packmol_input(system::MembraneProtein; input = "membrane_protein.inp", output = "membrane_protein.pdb")

Writes the input file of Packmol for a membrane protein system.
"""
function write_packmol_input(
    system::MembraneProtein;
    input = "membrane_protein.inp",
    output = "membrane_protein.pdb",
    debug = false,
)
    input, output = expanduser(input), expanduser(output)
    setup = _setup(system)
    (; lipid_placements, solvent_placements, lo, hi, summary) = setup
    println(summary)
    open(input, "w") do io
        println(io, "#\n# Packmol input file\n#\n# Generated by Packmol.jl\n#")
        foreach(line -> println(io, "# $line"), split(summary, "\n"))
        println(io, """
            #
            tolerance 2.0
            output $output
            filetype pdb
            seed -1
            pbc $(lo[1]) $(lo[2]) $(lo[3]) $(hi[1]) $(hi[2]) $(hi[3])

            structure $(system.protein_pdbfile)
                number 1
                fixed 0. 0. 0. 0. 0. 0.
            end structure
            """)
        foreach(lp -> _write_block(io, lp), lipid_placements)
        foreach(sp -> _write_block(io, sp), solvent_placements)
    end
    print(chomp("""
        Wrote file: $input

        ==================================================================
        """))
    return debug ? setup : nothing
end

@eval function packmol(system::MembraneProtein;
    output = "membrane_protein.pdb",
    $(map(_kwarg_signature, _packmol_recipe_kwargs)...),
)
    (; lipid_placements, solvent_placements, unitcell, center) = _setup(system)
    structure_types = StructureType{3,Float64}[
        structure_type(system.protein_pdbfile; number = 1, fixed = ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])),
    ]
    foreach(lp -> _push_structure_types!(structure_types, lp), lipid_placements)
    foreach(sp -> _push_structure_types!(structure_types, sp), solvent_placements)
    packmol_system = PackmolSystem(structure_types;
        output, tolerance = 2.0, add_box_sides = true, seed = -1,
        unitcell, unitcell_center = SVector{3,Float64}(center),
    )
    return packmol(packmol_system; $(_packmol_recipe_kwarg_names...))
end

"""
    packmol(system::MembraneProtein; output = "membrane_protein.pdb", kwargs...)

Builds and packs a membrane protein system, without writing an input file (as
[`write_packmol_input`](@ref write_packmol_input(::MembraneProtein)) followed by `packmol`).
Any other keyword is forwarded to the packing engine; returns the `PackmolSystem`.
"""
packmol(::MembraneProtein)
