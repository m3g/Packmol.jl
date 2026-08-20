@testitem "Molarity Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(1.0u"mol/L", "Molarity" => "Molarity") == 1.0u"mol/L"
    @test cconvert(1.0, "M" => "mol/L") == 1.0u"mol/L" # Number input

    # --- Molarity to Molality ---
    # Simple case: C=1, Ms=100, rho=1 -> m = 1/(1-0.1) = 1/0.9 = 1.111...
    @test cconvert(1.0u"mol/L", "M" => "mol/kg"; M_solute=100u"g/mol", rho_solution=1.0u"kg/L") ≈ 1.11u"mol/kg" rtol=1e-3
    # Ethanol case (~10 M -> ~22.76 m)
    @test cconvert(10.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ 26.019u"mol/kg" rtol=1e-3
    # Ethanol case (1 M -> ~1.036 m)
    @test cconvert(1.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 1.0684u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "mol/kg"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "mol/kg"; rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0u"mol/kg"

    # --- Molarity to MoleFraction ---
    # Ethanol case (1 M -> ~0.0183)
    @test cconvert(1.0u"M", "M" => "MoleFraction"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_1M) ≈ 0.01888 rtol=1e-3
    # Ethanol case (~10 M -> ~0.359)
    @test cconvert(10.0u"M", "M" => "MoleFraction"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_10M) ≈ 0.31914 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solute=M_EtOH, rho_solution=rho_sol_1M)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solvent=M_H2O, rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_1M) == 0.0

    # --- Molarity to MassFraction ---
    # Ethanol case (1 M -> ~0.04686)
    @test cconvert(1.0u"M", "M" => "MassFraction"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ (1.0*ustrip(M_EtOH)/1000)/ustrip(rho_sol_1M) # 0.046068 / 0.983
    @test cconvert(1.0u"M", "M" => "MassFraction"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 0.04691 rtol=1e-2
    # Ethanol case (~10 M -> ~0.5118)
    @test cconvert(10.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ (10.0*ustrip(M_EtOH)/1000)/ustrip(rho_sol_10M) # 0.46068 / 0.90
    @test cconvert(10.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ 0.54518 rtol=1e-2
    # Output as Percentage
    @test cconvert(1.0u"M", "M" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 4.69124 rtol=1e-2
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "w/w"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "w/w"; rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0
    @test cconvert(0.0u"M", "M" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0

    # --- Molarity to VolumeFraction ---
    # Ethanol case (1 M -> ~0.05839)
    @test cconvert(1.0u"M", "M" => "VolumeFraction"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ (1.0*ustrip(M_EtOH)/1000)/ustrip(rho_EtOH_pure) # 0.046068/0.789
    @test cconvert(1.0u"M", "M" => "VolumeFraction"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.05868 rtol=1e-3
    # Ethanol case (~10 M -> ~0.5839)
    @test cconvert(10.0u"M", "M" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.5868 rtol=1e-3
    # Output as Percentage
    @test cconvert(10.0u"M", "M" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 58.68 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "v/v"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "v/v"; rho_solute=rho_EtOH_pure)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0
    @test cconvert(0.0u"M", "M" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0

    # --- Molarity to NumberDensity ---
    @test cconvert(1.0u"mol/L", "M" => "Å^-3") ≈ 1.0u"mol/L" * Na rtol=1e-6
    @test cconvert(1.0u"mol/L", "M" => "Å^-3") ≈ 0.000602214u"Å^-3"  rtol=1e-6 
    # Missing Kwargs (none needed)
    @test cconvert(1.0u"M", "M" => "Å^-3") isa Quantity
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "Å^-3") == 0.0u"Å^-3"

end

@testitem "Molality Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(1.0u"mol/kg", "Molality" => "Molality") == 1.0u"mol/kg"
    @test cconvert(1.0, "mol/kg" => "molal") == 1.0u"mol/kg" # Number input

    # --- Molality to Molarity ---
    # Ethanol case (1 m -> ~0.97 M) - use density for 1m solution
    @test cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_1m) ≈ 0.9483u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "M"; rho_solution=rho_sol_1m)
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/L" # Use water density for zero conc

    # --- Molality to MoleFraction ---
    # Ethanol case (1 m -> ~0.0177)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) ≈ 1.0 / (1.0 + 1000/ustrip(M_H2O))
    @test cconvert(1.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) ≈ 0.017696 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "chi")
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) == 0.0

    # --- Molality to MassFraction ---
    # Ethanol case (1 m -> ~0.0440)
    m_s = 1.0 * ustrip(M_EtOH)
    mf = m_s / (1000 + m_s)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) ≈ mf
    @test cconvert(1.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) ≈ 0.044047 rtol=1e-3
    # Output Percentage
    @test cconvert(1.0u"mol/kg", "mol/kg" => "%m/m"; M_solute=M_EtOH) ≈ 4.4047 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "w/w")
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) == 0.0
    @test cconvert(0.0u"mol/kg", "mol/kg" => "%m/m"; M_solute=M_EtOH) == 0.0

    # --- Molality to VolumeFraction --- (Requires rho_solution, rho_solute)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) ≈ 0.0557 rtol=1e-3
    # Output Percentage
    @test cconvert(1.0u"mol/kg", "mol/kg" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) ≈ 5.5652 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solution=rho_sol_1m) # Missing rho_solute
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) # Missing M_solute
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0

    # --- Molality to NumberDensity --- (Requires rho_solution)
    # Ethanol case (1m -> intermediate molarity -> number density)
    nd = cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_sol_1m)
    C_interim = cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_1m)
    @test nd ≈ cconvert(C_interim, "mol/L" => "Å^-3")
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; rho_solution=rho_sol_1m) # Missing M_solute
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "MoleFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "MoleFraction" => "MoleFraction") == 0.5
    @test cconvert(0.5, "chi" => "x") == 0.5 # String variations

    # --- MoleFraction to Molarity ---
    # Ethanol case (χ=0.2 -> ~5.71 M) - use density for χ=0.2 solution
    @test cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ 7.88974u"mol/L" rtol=1e-3
    @test cconvert(1.0, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_EtOH_pure) ≈ 17.040u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_chi02)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solvent=M_H2O, rho_solution=rho_sol_chi02)
    # Edge Cases
    @test cconvert(0.0, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0u"mol/L"

    # --- MoleFraction to Molality ---
    # Ethanol case (χ=0.2 -> ~13.88 m)
    @test cconvert(0.2, "chi" => "mol/kg"; M_solvent=M_H2O) ≈ (0.2 / (0.8 * ustrip(M_H2O)/1000.0)) * u"mol/kg"
    @test cconvert(0.2, "chi" => "mol/kg"; M_solvent=M_H2O) ≈ 13.877u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "mol/kg")
    # Edge Cases
    @test cconvert(0.0, "chi" => "mol/kg"; M_solvent=M_H2O) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "chi" => "mol/kg"; M_solvent=M_H2O) # Requires chi < 1

    # --- MoleFraction to MassFraction ---
    # Ethanol case (χ=0.2 -> ~0.390)
    m_etoh = 0.2 * ustrip(M_EtOH)
    m_h2o = 0.8 * ustrip(M_H2O)
    mf = m_etoh / (m_etoh + m_h2o)
    @test cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ mf
    @test cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.38998 rtol=1e-3
    # Output Percentage
    @test cconvert(0.2, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 38.998 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "w/w"; M_solvent=M_H2O)
    # Edge Cases
    @test cconvert(0.0, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) == 0.0
    @test cconvert(1.0, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0
    @test cconvert(1.0, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O) == 100.0

    # --- MoleFraction to VolumeFraction --- (Requires densities)
    # Ethanol case (χ=0.2 -> ~0.265) - Requires intermediate Molarity calc
    @test cconvert(0.2, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ 0.4630 rtol=1e-3
    # Output Percentage
    @test cconvert(0.2, "chi" => "%v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ 46.30 rtol=1e-3 
    # Missing Kwargs (will cascade from intermediate call)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge cases
    @test cconvert(0.0, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    # Pure solute vf=1 needs density of pure solute as rho_solution
    @test cconvert(1.0, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) ≈ 1.0

     # --- MoleFraction to NumberDensity --- (Requires densities)
     # Ethanol case (χ=0.2 -> ~3.44e-3 Å^-3) - Requires intermediate Molarity calc
     @test cconvert(0.2, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ 0.0047513u"Å^-3" rtol=1e-3
     # Missing Kwargs
     @test_throws UndefKeywordError cconvert(0.2, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O) # Missing rho_solution
     # Edge Cases
     @test cconvert(0.0, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "MassFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "MassFraction" => "MassFraction") == 0.5
    @test cconvert(50.0, "%m/m" => "w/w") == 0.5 # Percent input, fraction output
    @test cconvert(0.5, "w/w" => "%m/m") == 50.0 # Fraction input, percent output

    # --- MassFraction to Molarity ---
    # Ethanol case (50% w/w -> ~9.92 M)
    @test cconvert(50.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3
    @test cconvert(0.5, "w/w" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "M"; rho_solution=rho_sol_50ww)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "M"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/L"
    @test cconvert(100.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_EtOH_pure) ≈ uconvert(u"mol/L", rho_EtOH_pure / M_EtOH) # Molarity of pure ethanol
    @test cconvert(100.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_EtOH_pure) ≈ 17.040u"mol/L" rtol=1e-3

    # --- MassFraction to Molality ---
    # Ethanol case (50% w/w -> 21.71 m)
    @test cconvert(50.0, "%m/m" => "mol/kg"; M_solute=M_EtOH) ≈ (0.5 / ustrip(M_EtOH)) / 0.5 * 1000 * u"mol/kg" # (m_s/M_s)/m_sv
    @test cconvert(50.0, "%m/m" => "mol/kg"; M_solute=M_EtOH) ≈ 21.707u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "mol/kg")
    # Edge Cases
    @test cconvert(0.0, "w/w" => "mol/kg"; M_solute=M_EtOH) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "w/w" => "mol/kg"; M_solute=M_EtOH) # mf must be < 1
    @test_throws ArgumentError cconvert(100.0, "%m/m" => "mol/kg"; M_solute=M_EtOH)

    # --- MassFraction to MoleFraction ---
    # Ethanol case (50% w/w -> χ ~ 0.259)
    n_etoh = 50 / ustrip(M_EtOH)
    n_h2o = 50 / ustrip(M_H2O)
    chi = n_etoh / (n_etoh + n_h2o)
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ chi
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.28111 rtol=1e-3
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.28111 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "chi"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "chi"; M_solvent=M_H2O)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 0.0
    @test cconvert(1.0, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0
    @test cconvert(100.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0

    # --- MassFraction to VolumeFraction ---
    # Ethanol case (50% w/w -> vf ~ 0.584)
    vf = (0.5 / ustrip(rho_EtOH_pure)) / (1.0 / ustrip(rho_sol_50ww))
    @test cconvert(50.0, "%m/m" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ vf
    @test cconvert(50.0, "%m/m" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 0.5822 rtol=1e-3 
    @test cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ (0.5/ustrip(rho_EtOH_pure))/(1.0/ustrip(rho_sol_50ww))
    @test cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 0.5822 rtol=1e-3
    # Output Percentage
    @test cconvert(50.0, "%m/m" => "%v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 58.216 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "v/v"; rho_solution=rho_sol_50ww)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) == 1.0

    # --- MassFraction to NumberDensity ---
    # Ethanol case (50% w/w -> ~5.97e-3 Å^-3) - intermediate molarity
    @test cconvert(50.0, "%m/m" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 0.00597u"Å^-3"  rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Å^-3"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Å^-3"; rho_solution=rho_sol_50ww) # Missing M_solute
    # Edge Cases
    @test cconvert(0.0, "w/w" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "VolumeFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "VolumeFraction" => "VolumeFraction") == 0.5
    @test cconvert(50.0, "%v/v" => "v/v") == 0.5 # Percent input, fraction output
    @test cconvert(0.5, "v/v" => "%v/v") == 50.0 # Fraction input, percent output

    # --- VolumeFraction to Molarity ---
    # Ethanol case (vf=0.5 -> ~8.56 M)
    @test cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ (0.5 * ustrip(rho_EtOH_pure)) / (ustrip(M_EtOH)/1000.0) * u"mol/L"
    @test cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 8.520u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "M"; rho_solute=rho_EtOH_pure)
    # Edge Cases
    @test cconvert(0.0, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0u"mol/L"
    @test cconvert(1.0, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ uconvert(u"mol/L", rho_EtOH_pure / M_EtOH) # Molarity of pure

    # --- VolumeFraction to MassFraction ---
    # Ethanol case (vf=0.5 -> mf ~ 0.424) - use rho for vf=0.5
    mf = (0.5 * ustrip(rho_EtOH_pure)) / ustrip(rho_sol_50vv)
    @test cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ mf
    @test cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 0.4220 rtol=1e-3
    # Output Percentage
    @test cconvert(50.0, "%v/v" => "%m/m"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 42.20 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "w/w"; rho_solution=rho_sol_50vv)
    # Edge Cases
    @test cconvert(0.0, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) == 1.0

    # --- VolumeFraction to Molality --- (Requires intermediate Molarity/MassFraction)
    # Ethanol case (vf=0.5 -> ~15.6 m)
    @test cconvert(0.5, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 15.85u"mol/kg" rtol=1e-3 # Chain M->m or MF->m
    # Missing Kwargs (will cascade)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) # Fails at MF->m step (mf=1)

    # --- VolumeFraction to MoleFraction --- (Requires intermediate Molarity/MassFraction)
    # Ethanol case (vf=0.5 -> χ ~ 0.218)
    @test cconvert(0.5, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 0.2221 rtol=1e-3
    # Missing Kwargs (will cascade)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0 - eps(1.0), "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) ≈ 1.0

    # --- VolumeFraction to NumberDensity ---
    # Ethanol case (vf=0.5 -> ~5.16e-3 Å^-3) - Intermediate Molarity
    @test cconvert(0.5, "v/v" => "Å^-3"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.0051308u"Å^-3"  rtol=1e-4
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "Å^-3"; M_solute=M_EtOH) # rho_solute missing
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "Å^-3"; rho_solute=rho_EtOH_pure) # M_solute missing
    # Edge Cases
    @test cconvert(0.0, "v/v" => "Å^-3"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0u"Å^-3"

end

@testitem "NumberDensity Conversions" setup=[CConvert] begin

    N_test = 0.005u"Å^-3" # Approx 8.3 M
    C_equiv = uconvert(u"mol/L", N_test / Na) # ~8.303 M
    # Need rho for this molarity, roughly rho_sol_10M=0.90 is closest defined
    rho_sol_Ntest = rho_sol_10M

    # --- Identity ---
    @test cconvert(N_test, "Å^-3" => "Å^-3") == N_test

    # --- NumberDensity to Molarity ---
    @test cconvert(N_test, "Å^-3" => "M") ≈ C_equiv  rtol=1e-6
    @test cconvert(N_test, "Å^-3" => "M") ≈ 8.303u"mol/L"  rtol=1e-4
    # Missing Kwargs (none needed)
    @test cconvert(N_test, "Å^-3" => "M") isa Quantity
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "M") == 0.0u"mol/L"

    # --- NumberDensity to Molality --- (Intermediate Molarity)
    # Ethanol case (N_test -> ~16.8 m)
    @test cconvert(N_test, "Å^-3" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 17.95u"mol/kg" rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "mol/kg"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "mol/kg"; rho_solution=rho_sol_Ntest) # Missing M_solute
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/kg"

    # --- NumberDensity to MoleFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> χ ~ 0.268)
    @test cconvert(N_test, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_Ntest) ≈ 0.2443 rtol=1e-3
     # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0

    # --- NumberDensity to MassFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> mf ~ 0.425)
    @test cconvert(N_test, "Å^-3" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 0.452649 rtol=1e-3
    # Output Percentage
    @test cconvert(N_test, "Å^-3" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 45.26 rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "w/w"; M_solute=M_EtOH) # Missing rho_solution
     # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "w/w"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0

    # --- NumberDensity to VolumeFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> vf ~ 0.485)
    @test cconvert(N_test, "Å^-3" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.4872 rtol=1e-3
    # Output Percentage
    @test cconvert(N_test, "Å^-3" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 48.72 rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "v/v"; M_solute=M_EtOH) # Missing rho_solute
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0

end

@testitem "Wrapper Functionality" setup=[CConvert] begin

    # Test % handling
    @test cconvert(50.0, "%m/m" => "w/w") == 0.5
    @test cconvert(0.5, "w/w" => "%m/m") == 50.0
    @test cconvert(50.0, "%m/m" => "MassFraction") == 0.5 # Target type name
    @test cconvert(0.5, "MassFraction" => "%m/m") == 50.0 # Source type name

    # Error fraction input > 1
    @test_throws ArgumentError cconvert(1.1, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test cconvert(50.0, " %m/m " => " Molarity "; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3

    # Test unknown units
    @test_throws ArgumentError cconvert(1.0, "bad_unit" => "Molarity")
    @test_throws ArgumentError cconvert(1.0, "Molarity" => "bad_unit")
    @test_throws ArgumentError cconvert(1.0u"mol/kg", "Molarity" => "Molarity")

    # Test missing kwargs via wrapper
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Molarity"; M_solute=M_EtOH) # Missing rho_solution

end

@testitem "Round Trips" setup=[CConvert] begin

    # Use χ=0.2 point
    chi1 = 0.2
    M1 = cconvert(chi1, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02)
    m1 = cconvert(chi1, "chi" => "mol/kg"; M_solvent=M_H2O)
    mf1_perc = cconvert(chi1, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O)
    vf1_perc = cconvert(chi1, "chi" => "%v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02)
    N1 = cconvert(chi1, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02)

    # Test round trips back to mole fraction (expect approx equality)
    @test cconvert(M1, "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3
    @test cconvert(m1, "mol/kg" => "chi"; M_solvent=M_H2O) ≈ chi1 rtol=1e-3
    @test cconvert(mf1_perc, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ chi1 rtol=1e-3
    @test cconvert(vf1_perc, "%v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3
    @test cconvert(N1, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3

end
