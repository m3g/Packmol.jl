@testitem "Membrane constructors" setup=[TestFiles] begin
    using Packmol
    using Unitful

    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    water = TestFiles.test_file("structure_files/water.pdb")

    # Built-in lipid: only the composition and the size
    sys = Membrane(lipids="POPC", total_lipids=20)
    @test sys.type == :bilayer
    @test sys.lipid_names == ["POPC"]
    @test basename(only(sys.lipid_pdbfiles)) == "popc.pdb"
    @test sys.lipid_head == [1]
    @test sys.lipid_tail == [[88, 131]]
    @test sys.area_per_lipid ≈ 1.05 * 64.0u"Å^2"
    @test sys.flexibility == 0.15
    @test basename(only(sys.solvent_pdbfiles)) == "water.pdb"
    @test sys.solvent_density == 1.0u"g/mL"
    @test sys.solvent_layer_width == 15.0u"Å"
    @test 25.0u"Å" < only(sys.lipid_head_tail_length) < 27.0u"Å"
    # case insensitive, mixture: area from the molar average
    sys = Membrane(lipids=["popc", "Pope"], lipid_molar_ratio=[3, 1], total_lipids=20)
    @test sys.lipid_names == ["POPC", "POPE"]
    @test sys.lipid_weight ≈ [0.75, 0.25]
    @test sys.area_per_lipid ≈ 1.05 * (0.75 * 64.0 + 0.25 * 57.0) * u"Å^2"
    # equal amounts by default
    @test Membrane(lipids=["DPPC", "DOPC"], total_lipids=20).lipid_weight ≈ [0.5, 0.5]

    # Lipid given as a file: head, chain ends (one or several) and area
    sys = Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0u"g/mL",
    )
    @test sys.lipid_pdbfiles == [lipid]
    @test sys.solvent_pdbfiles == [water]
    @test sys.total_lipids == 20
    @test isnothing(sys.total_area)
    @test sys.lipid_tail == [[1]]
    @test sys.lipid_head_tail_length[1] ≈ 11.574471473693848u"Å" rtol = 1e-6
    sys2 = Membrane(lipids=lipid, lipid_head=[31], lipid_tail=[[1, 2]], area_per_lipid=60.0, total_lipids=20)
    @test sys2.lipid_tail == [[1, 2]]
    @test sys2.lipid_head_tail_length[1] < sys.lipid_head_tail_length[1] + 2u"Å"
    # a file and a built-in lipid together
    sys = Membrane(
        lipids=["POPC", lipid], lipid_head=[1, 31], lipid_tail=[[88, 131], 1], area_per_lipid=60.0,
        total_lipids=20,
    )
    @test sys.lipid_names == ["POPC", "palmitoil.pdb"]

    # neither a built-in lipid nor a file
    @test_throws ArgumentError Membrane(lipids="XYZC", total_lipids=20)
    # a file without its head and chain ends, or without the area
    @test_throws ArgumentError Membrane(lipids=lipid, area_per_lipid=60.0, total_lipids=20)
    @test_throws ArgumentError Membrane(lipids=lipid, lipid_head=[31], lipid_tail=[1], total_lipids=20)
    @test_throws ArgumentError Membrane(lipids=lipid, lipid_head=[31], area_per_lipid=60.0, total_lipids=20)
    # invalid type
    @test_throws ArgumentError Membrane(type=:trilayer, lipids="POPC", total_lipids=20)
    # mismatched lengths
    @test_throws ArgumentError Membrane(lipids=["POPC", "POPE"], lipid_molar_ratio=[1.0], total_lipids=20)
    @test_throws ArgumentError Membrane(
        lipids=[lipid, lipid], lipid_head=[31], lipid_tail=[1], area_per_lipid=60.0, total_lipids=20,
    )
    # neither / both of total_area, total_lipids given
    @test_throws ArgumentError Membrane(lipids="POPC")
    @test_throws ArgumentError Membrane(lipids="POPC", total_lipids=20, total_area=2000.0)
    # out-of-range / degenerate head and chain-end indices
    @test_throws ArgumentError Membrane(lipids=lipid, lipid_head=[9999], lipid_tail=[1], area_per_lipid=60.0, total_lipids=20)
    @test_throws ArgumentError Membrane(lipids=lipid, lipid_head=[31], lipid_tail=[[1, 9999]], area_per_lipid=60.0, total_lipids=20)
    @test_throws ArgumentError Membrane(lipids=lipid, lipid_head=[1], lipid_tail=[1], area_per_lipid=60.0, total_lipids=20)
    # flexibility out of (0, 1)
    @test_throws ArgumentError Membrane(lipids="POPC", total_lipids=20, flexibility=1.0)
    @test_throws ArgumentError Membrane(lipids="POPC", total_lipids=20, flexibility=0.0)

    # density without units warns and assumes g/mL
    sys2 = @test_logs (:warn,) Membrane(lipids="POPC", total_lipids=20, solvent=water, solvent_density=1.0)
    @test sys2.solvent_density == 1.0u"g/mL"
end

@testitem "Membrane show method" setup=[TestFiles] begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    water = TestFiles.test_file("structure_files/water.pdb")
    sys = Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0u"g/mL",
    )
    @test parse_show(sys) ≈ """
    ==================================================================
    Membrane properties (bilayer):
    ==================================================================
        Lipids: palmitoil.pdb
        Lipid molar ratio (normalized): [1.0]
        Lipid head-to-tail lengths: 11.57 Å
        Solvents: water.pdb
        Solvent molar ratio (normalized): [1.0]
        Solvent density: 1.0 g mL^-1
        Solvent layer width: 15.0 Å
        Area per lipid: 60.0 Å^2
        Total lipids: 20
        Flexibility: 0.15
    ==================================================================
    """
end

@testitem "Membrane bilayer geometry" begin
    using Packmol
    using Unitful

    sys = Membrane(lipids=["POPC", "DOPC"], total_lipids=20)
    tmp_input_file = tempname() * ".inp"
    lipid_placements, solvent_placements, unitcell = redirect_stdout(devnull) do
        write_packmol_input(sys; input=tmp_input_file, debug=true)
    end
    @test isfile(tmp_input_file)

    L = ustrip.(u"Å", sys.lipid_head_tail_length)
    head_planes = [6.0 + 0.85 * l for l in L]
    zsolvent = maximum(head_planes) + 2.0
    a, b, c, α, β, γ = Packmol._unitcell_abc_angles(unitcell)
    @test a ≈ b
    @test a * b ≈ 10 * ustrip(u"Å^2", sys.area_per_lipid)
    @test c ≈ 2 * (zsolvent + 15.0)
    @test α ≈ β ≈ γ ≈ 90.0

    @test length(lipid_placements) == 4 # 2 lipid types x 2 leaflets
    @test length(solvent_placements) == 2 # 1 solvent x 2 slabs
    @test sum(lp.number for lp in lipid_placements) == 20
    upper = filter(lp -> lp.upper, lipid_placements)
    lower = filter(lp -> !lp.upper, lipid_placements)
    @test [lp.head_plane for lp in upper] ≈ head_planes
    @test [lp.head_plane for lp in lower] ≈ -head_planes
    @test all(lp.tail_plane == 6.0 for lp in upper)
    @test all(lp.tail_plane == -6.0 for lp in lower)
    @test all(isnothing(lp.floor) for lp in lipid_placements)
    @test upper[1].tail_idx == [88, 131]
    @test solvent_placements[1].zlo ≈ zsolvent
    @test solvent_placements[1].zhi ≈ c / 2
    @test solvent_placements[2].zhi ≈ -zsolvent

    # The generated input file is valid for the native engine
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == sum(lp.number for lp in lipid_placements) + sum(sp.number for sp in solvent_placements)
    rm(tmp_input_file, force=true)
end

@testitem "Membrane monolayer geometry" begin
    using Packmol
    using Unitful

    sys = Membrane(type=:monolayer, lipids="DPPC", total_lipids=10)
    tmp_input_file = tempname() * ".inp"
    lipid_placements, solvent_placements, unitcell = redirect_stdout(devnull) do
        write_packmol_input(sys; input=tmp_input_file, debug=true)
    end
    rm(tmp_input_file, force=true)

    a, b, c, = Packmol._unitcell_abc_angles(unitcell)
    L = ustrip(u"Å", only(sys.lipid_head_tail_length))
    zsolvent = 6.0 + 0.85 * L + 2.0
    @test c ≈ zsolvent + 15.0 # 1 leaflet + 1 solvent layer
    leaflet = only(lipid_placements)
    slab = only(solvent_placements)
    @test leaflet.upper
    @test leaflet.floor ≈ -c / 2
    @test leaflet.tail_plane ≈ -c / 2 + 6.0
    @test slab.zlo ≈ -c / 2 + zsolvent
    @test slab.zhi ≈ c / 2
end

@testitem "Membrane multi-lipid, multi-solvent, total_area" setup=[TestFiles] begin
    using Packmol
    using Unitful

    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    water = TestFiles.test_file("structure_files/water.pdb")
    diatomic = TestFiles.test_file("structure_files/diatomic.pdb")

    sys = Membrane(
        lipids=[lipid, lipid], lipid_head=[31, 30], lipid_tail=[1, 2], lipid_molar_ratio=[2.0, 1.0],
        area_per_lipid=60.0, total_area=2000.0,
        solvent=[water, diatomic], solvent_molar_ratio=[10.0, 1.0],
        solvent_layer_width=12.0, solvent_density=1.0u"g/mL",
    )
    @test sys.lipid_weight ≈ [2 / 3, 1 / 3]
    @test sys.solvent_weight ≈ [10 / 11, 1 / 11]

    tmp_input_file = tempname() * ".inp"
    lipid_placements, solvent_placements, unitcell = redirect_stdout(devnull) do
        write_packmol_input(sys; input=tmp_input_file, debug=true)
    end
    rm(tmp_input_file, force=true)

    a, b, c, = Packmol._unitcell_abc_angles(unitcell)
    @test a * b ≈ 2000.0 rtol = 1e-6

    @test length(lipid_placements) == 4 # 2 lipid types x 2 leaflets
    @test length(solvent_placements) == 4 # 2 solvents x 2 slabs
    # counts sum exactly to the requested lipids_per_leaflet, per leaflet
    lipids_per_leaflet = round(Int, 2000.0 / 60.0)
    for leaflet_placements in (lipid_placements[1:2], lipid_placements[3:4])
        @test sum(lp.number for lp in leaflet_placements) == lipids_per_leaflet
    end
end

@testitem "Membrane packmol (bilayer)" begin
    using Packmol
    using PDBToolsBase: read_pdb
    using Unitful

    sys = Membrane(lipids="POPC", total_lipids=20)
    tmp_output_file = tempname() * ".pdb"
    packed_sys = packmol(sys; output=tmp_output_file, seed=42, iprint=0)
    @test packed_sys isa Packmol.PackmolSystem
    @test packed_sys.status == :packing_ready
    @test isfile(tmp_output_file)

    # The heads above their planes and the chain ends within 6 Å of the middle, in both
    # leaflets (the upper one first)
    atoms = read_pdb(tmp_output_file)
    natoms = 134
    head_plane = 6.0 + 0.85 * ustrip(u"Å", only(sys.lipid_head_tail_length))
    for i in 1:20
        mol = atoms[(i-1)*natoms+1:i*natoms]
        s = i <= 10 ? 1 : -1
        @test s * mol[1].z > head_plane - 0.5
        @test s * mol[88].z < 6.5
        @test s * mol[131].z < 6.5
    end
    rm(tmp_output_file, force=true)
end

@testitem "Membrane packmol (monolayer)" begin
    using Packmol

    sys = Membrane(type=:monolayer, lipids="DMPC", total_lipids=10)
    tmp_output_file = tempname() * ".pdb"
    packed_sys = packmol(sys; output=tmp_output_file, seed=42, iprint=0)
    @test packed_sys isa Packmol.PackmolSystem
    @test packed_sys.status == :packing_ready
    @test isfile(tmp_output_file)
    rm(tmp_output_file, force=true)
end
