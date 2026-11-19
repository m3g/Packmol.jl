@testitem "rotatable bonds: fragments, links, exclusions" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: read_packmol_input, _nmonomers, _molecule_bonds, _fragments
    # A chain of 30 beads, 1.54 Å apart (no CONECT records)
    file = PolymerFiles.polymer_file("rot_chain.pdb"; nmon = 10)
    input = joinpath(PolymerFiles.dir, "rot.inp")
    block(lines) = (write(input, "tolerance 2.0\noutput rot_out.pdb\nstructure rot_chain.pdb\n number 4\n$lines\nend structure\n");
        input)
    sys = read_packmol_input(block("bond_distance 1.6\nrotatable 10 11 20 21\ninside sphere 0. 0. 0. 30."))
    p = only(sys.polymers)
    @test _nmonomers(p) == 3
    @test p.monomer_atoms == [collect(1:10), collect(11:20), collect(21:30)]
    @test length(sys.structure_types) == 3
    @test all(st -> st.radscale == 1.0 && length(st.constraints) == 1, sys.structure_types)
    # each rotatable bond: the bond and the bonds to the neighbors of its ends
    @test length(only(p.links)) == 2 * 3
    @test all(l -> isapprox(l.d0, 1.54; atol = 2e-3) || isapprox(l.d0, 2.55; atol = 0.02), only(p.links))
    # the pairs of different fragments up to three bonds apart
    @test length(only(p.exclusions)) == 2 * 6
    @test (Int32(10), Int32(13)) in only(p.exclusions)
    @test !((Int32(10), Int32(14)) in only(p.exclusions))
    # errors
    @test_throws ArgumentError read_packmol_input(block("rotatable 10 11"))                     # no bonds
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrotatable 10 12"))  # not bonded
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrotatable 10"))     # not a pair
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6"))                   # no rotatable bonds
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrotatable 10 11\nrepeat 10\n monomer 3\n link 3 1\nend repeat"))
    # A branched molecule: POPC, with the bonds from the distances (no CONECT records)
    popc = joinpath(Packmol.RecipesDirectory, "lipids", "popc.pdb")
    st = structure_type(popc; number = 1)
    bonds = _molecule_bonds(st, 1.8)
    @test sum(length, bonds) ÷ 2 == 133 # a tree of 134 atoms (no pairs of hydrogens)
    _, fragments = _fragments(bonds, [(25, 28), (33, 45), (57, 60), (42, 92), (101, 104)])
    @test length.(fragments) == [27, 17, 15, 32, 12, 31]
    # a bond in a ring cannot rotate: a ring of the beads of the chain
    ring = joinpath(PolymerFiles.dir, "ring.pdb")
    open(ring, "w") do io
        for i in 1:6
            θ = 2π * i / 6
            println(io, rpad("ATOM  " * lpad(i, 5) * "  C   RNG A   1", 30),
                join(lpad(string(round(v, digits = 3)), 8) for v in (1.54 * cos(θ), 1.54 * sin(θ), 0.0)), "  1.00  0.00")
        end
    end
    @test_throws ArgumentError PackmolSystem([flexible_structure(ring; number = 1, rotatable = [(1, 2)], bond_distance = 1.6)];
        output = joinpath(PolymerFiles.dir, "ring_out.pdb"))
end

@testitem "packing of molecules with rotatable bonds" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: PDBToolsBase
    file = PolymerFiles.polymer_file("rot_chains.pdb"; nmon = 10, nconf = 2, seed = 3)
    # input file
    input = joinpath(PolymerFiles.dir, "rot_pack.inp")
    write(input, """
        tolerance 2.0
        output rot_pack_out.pdb
        pbc 30. 30. 30.
        structure rot_chains.pdb
          number 12
          bond_distance 1.6
          rotatable 10 11 20 21
        end structure
        """)
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    atoms = PDBToolsBase.read_pdb(joinpath(PolymerFiles.dir, "rot_pack_out.pdb")).atoms
    @test length(atoms) == 12 * 30
    @test length(unique(a.chain for a in atoms)) == 12
    bonds, _ = PolymerFiles.check_chains(atoms, 30)
    @test all(isapprox.(bonds, 1.54; atol = 0.1))
    # Julia API
    chains = flexible_structure(file; number = 12, rotatable = [(10, 11), (20, 21)], bond_distance = 1.6)
    sys = PackmolSystem([chains]; output = joinpath(PolymerFiles.dir, "rot_api_out.pdb"), tolerance = 2.0,
        unitcell = [30.0 0 0; 0 30.0 0; 0 0 30.0])
    @test length(sys.structure_types) == 2 * 3
    sys = redirect_stdout(devnull) do
        packmol(sys)
    end
    @test sys.status == :packing_ready
    @test_throws ArgumentError flexible_structure(file; number = 1, rotatable = Tuple{Int,Int}[])
end

@testitem "rotatable bonds in the monomers of polymers" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: read_packmol_input
    # A chain of 30 beads, 1.54 Å apart, in 10 monomers of 3 beads
    file = PolymerFiles.polymer_file("rot_polymer.pdb"; nmon = 10)
    input = joinpath(PolymerFiles.dir, "rot_polymer.inp")
    block(lines) = (write(input, "tolerance 2.0\noutput rot_polymer_out.pdb\nstructure rot_polymer.pdb\n number 4\n$lines\nend structure\n");
        input)
    sys = read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 3 1\nend repeat"))
    p = only(sys.polymers)
    # the same links and exclusions as the molecule cut at the same bonds
    bonds = join(("$(3m) $(3m+1)" for m in 1:9), " ")
    q = only(read_packmol_input(block("bond_distance 1.6\nrotatable $bonds")).polymers)
    @test p.monomer_atoms == q.monomer_atoms
    pairs(p) = Set((p.monomer_atoms[l.fa][l.a], p.monomer_atoms[l.fb][l.b], round(l.d0; digits = 6)) for l in only(p.links))
    @test length(only(p.links)) == 9 * 3
    @test pairs(p) == pairs(q)
    @test only(p.exclusions) == only(q.exclusions)
    # with links too, and in the Julia API
    sys = read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 3 1\n  link 2 2 3.0\nend repeat"))
    @test length(only(only(sys.polymers).links)) == 9 * 4
    chains = polymer_structure(file; number = 4, bond_distance = 1.6, sequence = [10 => [monomer(3; rotatable = [(3, 1)])]])
    s = PackmolSystem([chains]; output = joinpath(PolymerFiles.dir, "x.pdb"))
    @test pairs(only(s.polymers)) == pairs(p)
    # errors
    @test_throws ArgumentError read_packmol_input(block("repeat 10\n monomer 3\n  rotatable 3 1\nend repeat"))                    # no bonds
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 3 2\nend repeat"))  # not bonded
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 4 1\nend repeat"))  # not in the monomer
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 3 4\nend repeat"))  # not in the next
    @test_throws ArgumentError read_packmol_input(block("bond_distance 1.6\nrepeat 10\n monomer 3\n  rotatable 3\nend repeat"))    # not a pair
    @test_throws ArgumentError monomer(3; rotatable = [(4, 1)])
    # packing: the chains whole, the bonds kept
    write(input, "tolerance 2.0\noutput rot_polymer_out.pdb\nstructure rot_polymer.pdb\n number 4\n inside sphere 0. 0. 0. 25.\n" *
        " bond_distance 1.6\n repeat 10\n  monomer 3\n   rotatable 3 1\n end repeat\nend structure\n")
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    atoms = Packmol.PDBToolsBase.read_pdb(joinpath(PolymerFiles.dir, "rot_polymer_out.pdb")).atoms
    bonds, _ = PolymerFiles.check_chains(atoms, 30)
    @test all(isapprox.(bonds, 1.54; atol = 0.1))
end
