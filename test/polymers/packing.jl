@testitem "packing of polymers" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: PDBToolsBase
    file = PolymerFiles.polymer_file("chains.pdb"; nmon = 10, nconf = 3, seed = 2)
    dir = PolymerFiles.dir
    # A single chain: a solution from the start
    input = joinpath(dir, "single.inp")
    write(input, "tolerance 2.0\noutput single_out.pdb\n" *
        PolymerFiles.polymer_block(file; number = 1, extra = "inside sphere 0. 0. 0. 40."))
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    # Chains in a periodic cell: whole in the output, each with its own chain identifier
    input = joinpath(dir, "pbc.inp")
    write(input, "tolerance 2.0\noutput pbc_out.pdb\npbc 30. 30. 30.\nseed 3\n" *
        PolymerFiles.polymer_block(file; number = 12))
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    atoms = PDBToolsBase.read_pdb(joinpath(dir, "pbc_out.pdb")).atoms
    @test length(atoms) == 12 * 30
    @test length(unique(a.chain for a in atoms)) == 12
    @test [a.resnum for a in atoms[1:30]] == [(i - 1) ÷ 3 + 1 for i in 1:30]
    bonds, _ = PolymerFiles.check_chains(atoms, 30)
    @test all(isapprox.(bonds, 1.54; atol = 0.1))
    # the closest pair beyond the bonds, with the minimum image
    x = [[a.x, a.y, a.z] for a in atoms]
    image(v) = v .- 30 .* round.(v ./ 30)
    dmin = minimum(sqrt(sum(abs2, image(x[i] - x[j])))
        for i in eachindex(x) for j in i+1:length(x) if !((i - 1) ÷ 30 == (j - 1) ÷ 30 && j == i + 1))
    @test dmin > 1.9
end
