mutable struct SolutionBoxUSC <: SolutionBox
    concentration_units::String
    solute_pdbfile::String
    solvent_pdbfile::String
    cossolvent_pdbfile::String
    density_table::Matrix{Float64}
    solute_molar_mass::Float64
    solvent_molar_mass::Float64
    cossolvent_molar_mass::Float64
end

"""
    SolutionBoxUSC(; 
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        density_table::Matrix{Float64},
        concentration_units = "x",
        solute_molar_mass = nothing, # optional
        solvent_molar_mass = nothing, # optional
        cossolvent_molar_mass = nothing, # optional
    )

Setup a system composed of a solute (U) a solvent (S) and a cossolvent (C). 

The concentration units of the density table can be provided explicitly and
are assumed by default to be the molar fraction, `x`, of the cossolvent.

The molar masses of the solute, solvent, and cossolvent can be provided manually. If
not, they will be computed from the atom types in the PDB file.

"""
function SolutionBoxUSC(;
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        density_table::Matrix{Float64},
        concentration_units::Union{Nothing,String} = nothing,
        solute_molar_mass::Union{Nothing,Real} = nothing,
        solvent_molar_mass::Union{Nothing,Real} = nothing,
        cossolvent_molar_mass::Union{Nothing,Real} = nothing,
    )
    if isnothing(concentration_units)
        concentration_units = "x"
        @warn "Concentration units not provided, assuming molar fraction." _file=nothing _line=nothing
    end
    if !(density_table[begin, 1] == 0.0)
        throw(ArgumentError("First line of density table must be the density of pure solvent, with cossolvent concentration equal to 0.0"))
    end
    if concentration_units in ("x", "vv", "mm") && !(density_table[end, 1] == 1.0)
        throw(ArgumentError("Last line of density table must be the density of pure cossolvent, with cossolvent concentration equal to 1.0"))
    end
    if isnothing(solute_molar_mass)
        solute_molar_mass = mass(read_pdb(solute_pdbfile))
    end
    if isnothing(solvent_molar_mass)
        solvent_molar_mass = mass(read_pdb(solvent_pdbfile))
    end
    if isnothing(cossolvent_molar_mass)
        cossolvent_molar_mass = mass(read_pdb(cossolvent_pdbfile))
    end
    system = SolutionBoxUSC(
        concentration_units,
        solute_pdbfile,
        solvent_pdbfile,
        cossolvent_pdbfile,
        density_table,
        solute_molar_mass,
        solvent_molar_mass,
        cossolvent_molar_mass,
    )
    if concentration_units == "mol/L"
        x = convert_concentration(system, density_table[end, 1], "mol/L" => "x")
        if !(x ≈ 1)
            throw(ArgumentError(chomp("""
            Conversion of concentration of last line into molar fraction gives: $x which is different from 1.0. 

            The last line of the density_table must contain the density of pure cossolvent.

            """)))
        end
    end
    return system
end

function Base.show(io::IO, ::MIME"text/plain", system::SolutionBoxUSC)
    print(io, chomp("""
    ==================================================================
    SolutionBoxUSC properties (Solute + Solvent + Cossolvent):
    ==================================================================
        Solute pdb file: $(basename(system.solute_pdbfile))
        Solvent pdb file: $(basename(system.solvent_pdbfile))
        Cossolvent pdb file: $(basename(system.cossolvent_pdbfile))
        Density of pure solvent: $(density_pure_solvent(system)) g/mL
        Density of pure cossolvent: $(density_pure_cossolvent(system)) g/mL
        Molar masses: 
            solute: $(system.solute_molar_mass) g/mol
            solvent: $(system.solvent_molar_mass) g/mol
            cossolvent: $(system.cossolvent_molar_mass) g/mol
        Concentration units: $(system.concentration_units) ($(unit_name(system.concentration_units)))
        Cocentration range: $(first(system.density_table[:, 1])) - $(last(system.density_table[:, 1]))
    ==================================================================
    """))
end

"""
    convert_concentration(
        system::SolutionBoxUSC,
        input_concentration, 
        units
    )

Convert concentration from one unit to another. The input
concentration is given in `input_concentration`, and the unit conversion 
is given by `units` keyword, that can be one of the following pairs:

The supported concentration units are:

- `"mol/L"`: molarity
- `"x"`: molar fraction
- `"vv"`: volume fraction
- `"mm"`: mass fraction

Conversion among types consists in passing the `units` keyword argument,
which is a pair of the form `"from" => "to"`, where `"from"` and `"to"`
are one of the supported units.

# Example

For example, to convert from molarity to molar fraction, use:

```julia
convert_concentration(system, 55.5, "mol/L" => "x")
```
where `system` is a `SolutionBoxUSC` object, and `55.5` is the molarity.

"""
function convert_concentration(
    system::SolutionBoxUSC,
    input_concentration::Real, 
    units::Pair{String,String}
)

    (; solvent_molar_mass, cossolvent_molar_mass) = system

    # If the units didn't change, just return the input concentrations
    first(units) == last(units) && return input_concentration

    # Obtain density of the solution by interpolation
    if first(units) != system.concentration_units
        input_concentration_units = system.concentration_units
        convert_density_table!(system, first(units))
        ρ = interpolate_concentration(system, input_concentration)
        convert_density_table!(system, input_concentration_units)
    else
        ρ = interpolate_concentration(system, input_concentration)
    end

    ρc = density_pure_cossolvent(system) # density of the pure cossolvent
    ρw = density_pure_solvent(system) # density of the pure solvent
    Mw = solvent_molar_mass # molar mass of solvent 
    Mc = cossolvent_molar_mass # molar mass of the cossolvent

    # nc and nw are the molar concentrations
    if first(units) == "vv"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Volume fraction must be in the [0,1] range."))
        end
        # mL * g/ml / g/mol = mol 
        vv = input_concentration
        Nc = (ρc * vv) / Mc # number of cossolvent molecules in 1mL of ideal solution
        Nw = ρw * (1 - vv) / Mw # number of solvent molecules in 1mL of ideal solution
        if last(units) == "x"
            x = Nc / (Nc + Nw) # molar fraction
            return fixrange(x)
        end
        if last(units) == "mol/L"
            mt = Nc * Mc + Nw * Mw # mass of the solution (for 1 mol total)
            v = mt / (1000*ρ) # volume of the solution (for 1 mol total)
            return Nc / v # mol/L
        end
        if last(units) == "mm"
            return fixrange(Nc * Mc / (Nc * Mc + Nw * Mw))
        end
    end

    if first(units) == "x"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Molar fraction must be in the [0,1] range."))
        end
        x = input_concentration # molar fraction
        if last(units) == "mol/L"
            m = x * Mc + (1 - x) * Mw # g/mol: mass of the solution
            v = m / (ρ*1000) # L/mol: volume
            return x / v # mol/L: molarity of the solution 
        end
        if last(units) == "vv"
            vc = Mc *  x / ρc # Volume of x mols of pure cossolvent 
            vw = Mw * (1 - x) / ρw # Volume of (1-x) mols of pure solvent 
            vv = vc / (vc + vw) # volume fraction of cossolvent in ideal solution
            return fixrange(vv)
        end
        if last(units) == "mm"
            return fixrange(x * Mc / (x * Mc + (1 - x) * Mw))
        end
    end

    if first(units) == "mol/L"
        pure_c = 1000 * ρc / Mc  
        if !(0 <= input_concentration <= pure_c) && !(input_concentration ≈ pure_c)
            @show input_concentration
            throw(ArgumentError("Cossolvent molarity must be in the [0,$pure_c] range."))
        end
        nc = input_concentration / 1000
        nw = (ρ - nc * Mc) / Mw
        if last(units) == "x"
            return fixrange(nc / (nc + nw))
        end
        if last(units) == "vv"
            vc = nc * Mc / ρc
            vw = nw * Mw / ρw
            return fixrange(vc / (vc + vw))
        end
        if last(units) == "mm"
            return fixrange(nc * Mc / (nc * Mc + nw * Mw))
        end
    end

    if first(units) == "mm"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Mass fraction must be in the [0,1] range."))
        end
        mm = input_concentration # mass fraction
        Nc = mm * 1 / Mc # mol of ethanol in 1g
        Nw = (1 - mm) / Mw # mol of water in 1g
        if last(units) == "x"
            return fixrange(Nc / (Nc + Nw)) # molar fraction
        end
        if last(units) == "vv"
            vc = Nc * Mc  / ρc 
            vw = Nw * Mw / ρw
            return fixrange(vc / (vc + vw)) # volume fraction
        end
        if last(units) == "mol/L"
            v = 1 / ρ # volume of 1g of solution
            return 1000 * Nc / v # mol/L
        end
    end
end

"""
    convert_density_table!(system::SolutionBoxUSC, target_units)

Converts the density table of the system from one unit to another. Returns the 
input `system` with the density table converted to the new units.

The target units may be one of: `"mol/L"`, `"x"`, `"vv"`, `"mm"`.

## Example

```julia
convert_density_table!(system, "mol/L")
```

"""
function convert_density_table!(
    system::SolutionBoxUSC,
    target_units::String;
)
    current_units = system.concentration_units
    new_density_table = copy(system.density_table)
    for irow in eachindex(eachrow(new_density_table))
        cin = new_density_table[irow, 1]
        ρ = new_density_table[irow, 2]
        cout = convert_concentration(system, cin, current_units => target_units)
        new_density_table[irow, 1] = cout
    end
    system.density_table .= new_density_table
    system.concentration_units = target_units
    return system
end

"""
    write_packmol_input(
        system::SolutionBoxUSC;
        concentration::Real, 
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real,
        cubic::Bool = false,
    )

Function that generates an input file for Packmol. 

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternativelly, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

If `cubic` is set to true, the box will be cubic, and the box sides will be
equal in all 3 dimensions, respecting the minimum margin provided.

"""
function write_packmol_input(
    system::SolutionBoxUSC;
    concentration::Real, 
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Real},Nothing} = nothing,
    margin::Union{Real,Nothing} = nothing, 
    cubic::Bool = false,
    # testing option
    debug = false,
)

    (; solute_pdbfile,
       solvent_pdbfile,
       cossolvent_pdbfile,
       solute_molar_mass,
       solvent_molar_mass,
       cossolvent_molar_mass,
    ) = system

    # molar masses (g/mol)
    Mp = solute_molar_mass
    Mc = cossolvent_molar_mass
    Mw = solvent_molar_mass

    # Check consistency of the concentrations given
    cunit = system.concentration_units
    ρs = @view(system.density_table[:, 1])
    if cunit == "mol/L" && (first(ρs) < 1e-3 && last(ρs) ≈ 1)
        throw(ArgumentError("Concentrations in density table do not appear to be in mol/L."))
    end
    if (cunit in ("x", "mm", "vv")) && (any(x -> !(0 <= x <= 1), ρs)) 
        throw(ArgumentError("Concentrations in density table outside [0,1] range, and units are: $cunit"))
    end

    # Find the density corresponding to the target concentration
    ρ = interpolate_concentration(system, concentration)

    # Obtain the concentration in all units, for testing
    c_x = convert_concentration(system, concentration, cunit => "x")
    c_vv = convert_concentration(system, concentration, cunit => "vv")
    cc_mol = convert_concentration(system, concentration, cunit => "mol/L")
    cc_mm = convert_concentration(system, concentration, cunit => "mm")

    # aliases for clearer formulas
    ρc = density_pure_cossolvent(system) 
    ρw = density_pure_solvent(system) 

    # Convert cossolvent concentration in molecules/Å³
    cc = CMC * cc_mol

    # Set box side
    if isnothing(box_sides) && isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided."))
    elseif !isnothing(box_sides) && !isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided, but not both."))
    end
    solute_atoms = read_pdb(system.solute_pdbfile)
    solute_extrema = round.(maxmin(solute_atoms).xlength; digits=3)
    if !isnothing(margin)
        box_sides = solute_extrema .+ 2 .* margin
    end

    # Box volume (Å³)
    if !cubic
        vbox = box_sides[1] * box_sides[2] * box_sides[3]
    else
        max_side = maximum(box_sides)
        vbox = max_side^3
        box_sides .= max_side
    end

    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the solution
    vs = vbox - CMV * Mp / ρ

    # number of cossolvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, cc * vs)

    # Number of solvent molecules
    if nc != 0
        nw = round(Int, nc * (1 - c_x) / c_x)
    else
        nw = round(Int, vs * (ρw/CMV) / Mw)
    end

    # Final density of the solution (not inclusing solute volume)
    ρ = CMV * (Mc * nc + Mw * nw) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = 1000 * (nc / vs) * CMV

    # Final solvent concentration (mol/L)
    cw_f = 1000 * (nw / vs) * CMV

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $cc_mol mol/L
        (of cossolvent)      = $c_x molar fraction
                             = $c_vv volume fraction
                             = $cc_mm mass fraction 
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
        Periodic box = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 

        Density of solution = $ρ g/mL
        Solute molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol
        Solvent molar mass = $Mw g/mol

        Number of cossolvent ($(basename(cossolvent_pdbfile))) molecules = $nc 
        Number of solvent ($(basename(solvent_pdbfile))) molecules = $nw 

        Final cossolvent concentration = $cc_f mol/L
                                    = $(CMC*cc_f) molecules/Å³
        Final solvent concentration = $cw_f mol/L
                                    = $(CMC*cw_f) molecules/Å³
                                    
        Final volume fraction = $((nc * Mc / ρc)/((nc * Mc / ρc) + (nw * Mw / ρw)))
        Final molar fraction = $(nc/(nc+nw))

        Cubic box requested: $cubic

        ==================================================================
        """
    println(summary)

    open(input, "w") do io
        print(io,
            """
            # 
            # Packmol input file
            # 
            # Generated by MolSimToolkit.jl
            #
            """
        )
        for line in split(summary, "\n")
            println(io, "# $line")
        end
        println(io,
            """
            #
            tolerance 2.0
            output $output
            add_box_sides 1.0
            filetype pdb
            seed -1
            packall
            pbc $(join( -1.0*l, " ")), $(join( l, " "))

            structure $solute_pdbfile
                number 1
                center
                fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $solvent_pdbfile
                number $nw
            end structure
            """)
        if nc > 0
            println(io,
                """
                structure $cossolvent_pdbfile
                    number $nc
                end structure
                """)
        end
    end
    print(chomp(
        """
        Wrote file: $input

        ==================================================================
        """))
    
    if debug 
        return nw, nc, 2*l
    else
        return nothing
    end
end # function write_packmol_input
