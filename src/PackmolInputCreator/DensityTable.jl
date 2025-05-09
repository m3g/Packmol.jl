struct DensityTable
    concentration::Vector{typeof(1.0u"mol/L")}
    density::Vector{typeof(1.0u"g/mL")}
end

function DensityTable(
    concentration_units::String,
    concentration::AbstractVector{<:Number},
    density::AbstractVector{<:Number},
    M_solvent::Quantity,
    M_cossolvent::Quantity,
)
    dvec = typeof(1.0u"g/mL")[uconvert(u"g/mL", 1u"g/mL" * d) for d in density]
    cvec = typeof(1.0u"mol/L")[
        cconvert(c, concentration_units => "mol/L"; M_solvent, M_solute=M_cossolvent, rho_solution=d) 
        for (c,d) in zip(concentration,dvec)
    ]
    return DensityTable(cvec, dvec)
end

function Base.show(io::IO, ::MIME"text/plain", density_table::DensityTable)
    uc = unit(eltype(density_table.concentration))
    ud = unit(eltype(density_table.density))
    print(io, """
    ==================================================================
    Density table:
    ==================================================================
    Concentration ($uc) |      Density ($ud)
    """)
    for (c,d) in zip(density_table.concentration, density_table.density)
        @printf(io, "%15.3f%9s | %15.3f\n", ustrip(c), ' ', ustrip(d))
    end
    print(io, chomp("""
    ==================================================================
    """))
end

#
# Densities of pure solvent and cossolvent
#
function density_pure_solvent(density_table::DensityTable) 
    cfirst = first(density_table.concentration)
    if !(_round(cfirst) ≈ 0.0*cfirst)
        return missing # pure solvent density is not in the table
    end
    first(density_table.density)
end
density_pure_solvent(system::SolutionBox) = density_pure_solvent(system.density_table)

function density_pure_cossolvent(density_table::DensityTable, M_solvent::Quantity, M_solute::Quantity) 
    clast = last(density_table.concentration)
    clast = cconvert(clast, "mol/L" => "x"; 
        M_solvent = M_solvent,  
        M_solute = M_solute, 
        rho_solution = last(density_table.density),
    )
    if !(_round(clast) ≈ 1.0)
        return missing # pure cossolvent density is not in the table
    end
    last(density_table.density[end])
end
function density_pure_cossolvent(system::SolutionBox) 
    return density_pure_cossolvent(
        system.density_table,
        system.solvent_molar_mass,
        system.cossolvent_molar_mass,
    )
end

#
# Interpolate density
#
function interpolate_density(system::SolutionBox, concentration::Number, concentration_units::String)
    return interpolate_density(
        system.density_table, 
        system.solvent_molar_mass,
        system.cossolvent_molar_mass,
        concentration, 
        concentration_units,
    )
end

function interpolate_density(
    density_table::DensityTable,
    M_solvent::Quantity,
    M_solute::Quantity,
    concentration::Number,
    concentration_units::String,
)
    dtab = density_table.density
    ctab = density_table.concentration
    dpc = density_pure_cossolvent(density_table, M_solvent, M_solute)
    for (ic,c) in enumerate(ctab)
        cc = cconvert(c, "mol/L" => concentration_units; 
            M_solvent = M_solvent,
            M_solute = M_solute,
            rho_solution = dtab[ic],
            rho_solute = dpc,
        )
        if (ic == 1 && concentration < cc) || (ic+1 > length(ctab))
            cfirst = cconvert(ctab[1], "mol/L" => concentration_units; 
                M_solvent, 
                M_solute, 
                rho_solution = dtab[1],
                rho_solute = dpc,
            ) 
            clast = cconvert(ctab[end], "mol/L" => concentration_units; 
                M_solvent = M_solvent,
                M_solute = M_solute,
                rho_solution = dtab[end],
                rho_solute = dpc, 
            )
            throw(ArgumentError("""\n
                Concentration out of range of density table.
                In units of $(concentration_units): $(cfirst) - $(clast)
                Concentration provided: $(concentration)
            """))
        end
        cnext = cconvert(ctab[ic+1], "mol/L" => concentration_units; 
            M_solvent = M_solvent,
            M_solute = M_solute,
            rho_solution = dtab[ic+1],
            rho_solute = dpc, 
        )
        if cc <= concentration <= cnext
            return dtab[ic] + (dtab[ic+1] - dtab[ic]) * (concentration - cc) / (cnext - cc)
        end
    end
end
