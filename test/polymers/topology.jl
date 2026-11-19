@testitem "minimum image of long vectors" begin
    using Packmol
    using StaticArrays: SVector
    unitcell = [30.0 0.0 0.0; 0.0 30.0 0.0; 0.0 0.0 30.0]
    # vectors several cells long (between the monomers of an unwrapped chain)
    @test Packmol._min_image(SVector(-56.2, 0.5, 0.0), unitcell) ≈ SVector(3.8, 0.5, 0.0)
    @test Packmol._min_image(SVector(91.0, -62.0, 29.0), unitcell) ≈ SVector(1.0, -2.0, -1.0)
    triclinic = [30.0 10.0 0.0; 0.0 28.0 0.0; 0.0 0.0 25.0]
    v = SVector(1.0, 2.0, -1.5)
    @test Packmol._min_image(v + 3 * SVector(10.0, 28.0, 0.0) - 2 * SVector(30.0, 0.0, 0.0), triclinic) ≈ v
end

@testitem "polymers whole across several cells" setup=[PolymerFiles] begin
    using Packmol
    using Packmol: PDBToolsBase
    # A straight chain (about 37 Å) in a periodic cube of 15 Å: it crosses the faces at
    # least twice, and must be written whole
    file = PolymerFiles.polymer_file("straight.pdb"; nmon = 10, nconf = 1)
    dir = PolymerFiles.dir
    input = joinpath(dir, "straight.inp")
    write(input, "tolerance 2.0\noutput straight_out.pdb\npbc 15. 15. 15.\n" *
        PolymerFiles.polymer_block(file; number = 1))
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    atoms = PDBToolsBase.read_pdb(joinpath(dir, "straight_out.pdb")).atoms
    bonds, _ = PolymerFiles.check_chains(atoms, 30)
    @test all(isapprox.(bonds, 1.54; atol = 0.1))
end
