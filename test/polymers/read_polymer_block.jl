@testitem "polymer blocks of input files" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: read_packmol_input, _npolymers, _nmonomers, _nconformations, _monomer_type
    file = PolymerFiles.polymer_file("chain.pdb"; nmon = 10, nconf = 2)
    input = joinpath(PolymerFiles.dir, "chain.inp")
    write(input, """
        tolerance 2.0
        output chain_out.pdb
        $(PolymerFiles.polymer_block(file; number = 5, extra = "inside sphere 0. 0. 0. 30."))
        structure $(basename(file))
          number 2
          inside sphere 0. 0. 0. 30.
        end structure
        """)
    sys = read_packmol_input(input)
    @test length(sys.polymers) == 1
    p = only(sys.polymers)
    @test _nmonomers(p) == 10
    @test _nconformations(p) == 2
    @test _npolymers(p) == 5
    # a structure type for each monomer of each conformation, then the rigid molecules
    @test length(sys.structure_types) == 2 * 10 + 2
    @test sys.structure_types[_monomer_type(p, 2, 1)].number_of_molecules == 2
    @test all(st -> st.radscale == 1.0, sys.structure_types[1:20])
    @test all(st -> length(st.constraints) == 1, sys.structure_types[1:20])
    # the links of the last monomer of a linear polymer are ignored; d0 from the structure
    @test [(l.fa, l.fb) for l in p.links[1]] == [(m, m + 1) for m in 1:9]
    @test all(l -> isapprox(l.d0, 1.54; atol = 2e-3), p.links[1])
    # the bonded beads of consecutive monomers are excluded pairs
    @test (Int32(3), Int32(4)) in p.exclusions[1]
    @test !((Int32(2), Int32(4)) in p.exclusions[1])
    @test sys.nmols == 5 * 10 + 2
    # errors
    bad(block) = (write(input, "output x.pdb\nstructure $(basename(file))\n number 1\n$block\nend structure\n");
        @test_throws ArgumentError read_packmol_input(input))
    bad("repeat 10\n fragment 4\n link 4 1\nend repeat")         # 40 atoms, the structure has 30
    bad("repeat 10\n fragment 3\n link 3 4\nend repeat")         # atom 4 of the next monomer
    bad("repeat 10\n fragment 3\n link 4 1\nend repeat")         # atom 4 of this monomer
    bad("repeat 10\n fragment 3\n link 3 1")                     # not closed
    bad("link 3 1\n")                                           # not after a fragment line
    bad("repeat 10\n fragment 3\n link 3 1\nend repeat\nlink 3 1")  # after the repeat block
    bad("repeat 10\n fragment 3\nend repeat")                    # no links
    bad("cyclic\n")                                             # no fragments
    bad("bond_distance 1.6\nrotatable 3 4\nrepeat 10\n fragment 3\n link 3 1\nend repeat") # both
    # fragments outside repeat blocks: end groups, or a sequence without repeats
    write(input, "output x.pdb\nstructure $(basename(file))\n number 1\n fragment 3\n  link 3 1\n" *
        " repeat 8\n  fragment 3\n   link 3 1\n end repeat\n fragment 3\nend structure\n")
    q = only(read_packmol_input(input).polymers)
    @test _nmonomers(q) == 10
    @test [(l.fa, l.fb) for l in q.links[1]] == [(m, m + 1) for m in 1:9]
    write(input, "output x.pdb\nstructure $(basename(file))\n number 1\n" * "fragment 3\n link 3 1\n"^10 * "end structure\n")
    @test [(l.fa, l.fb) for l in only(read_packmol_input(input).polymers).links[1]] == [(m, m + 1) for m in 1:9]
    # cyclic: the last monomer linked to the first
    write(input, "output x.pdb\n" * PolymerFiles.polymer_block(file; number = 1, cyclic = true))
    @test [(l.fa, l.fb) for l in only(read_packmol_input(input).polymers).links[1]] == [[(m, m + 1) for m in 1:9]; (10, 1)]
end

@testitem "flexible_structure_type with fragments (Julia API)" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: _npolymers, _nmonomers
    file = PolymerFiles.polymer_file("chain_api.pdb"; nmon = 10, nconf = 2)
    chains = flexible_structure_type(file; number = 4, fragments = [10 => fragment(3; links = [(3, 1)])],
        constraints = [Packmol.InsideSphere([0.0, 0.0, 0.0], 30.0)])
    water = structure_type(file; number = 1)
    sys = PackmolSystem([chains, water]; output = joinpath(PolymerFiles.dir, "api.pdb"))
    p = only(sys.polymers)
    @test _nmonomers(p) == 10
    @test _npolymers(p) == 4
    @test sys.nmols == 4 * 10 + 1
    # links with the distance, constant and flat bottom given
    @test chains isa Packmol.FlexibleStructureType
    # single fragments and repeats, in any order
    ends = flexible_structure_type(file; number = 1,
        fragments = [fragment(3; links = [(3, 1)]), 8 => [fragment(3; links = [(3, 1)])], fragment(3)])
    @test length(ends.spec.fragments) == 10
    @test isempty(last(ends.spec.fragments).links)
    m = fragment(3; links = [(3, 1, 1.6, 5.0, 0.1)])
    @test only(m.links).d0 == 1.6 && only(m.links).k == 5.0 && only(m.links).w == 0.1
    @test_throws ArgumentError fragment(3; links = [(4, 1)])
    @test_throws ArgumentError fragment(3; links = [(1,)])
    @test_throws ArgumentError flexible_structure_type(file; number = 1, fragments = [9 => [fragment(3; links = [(3, 1)])]])
    @test_throws ArgumentError flexible_structure_type(file; number = 1)                        # neither
    @test_throws ArgumentError flexible_structure_type(file; number = 1, fragments = [10 => fragment(3)],
        rotatable = [(3, 4)])                                                                  # both
    @test_throws ArgumentError flexible_structure_type(file; number = 1, rotatable = [(3, 4)], cyclic = true)
    @test_throws ArgumentError flexible_structure_type(file; number = 1, fragments = [0 => fragment(3)])
    @test_throws ArgumentError flexible_structure_type(file; number = 1, fragments = [3])
end
