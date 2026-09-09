@testitem "Membrane constructors" begin
    using Packmol
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")

    # basic construction (single lipid, single solvent)
    sys = Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )
    @test sys.type == :bilayer # default
    @test sys.lipid_pdbfiles == [lipid]
    @test sys.solvent_pdbfiles == [water]
    @test sys.total_lipids == 20
    @test isnothing(sys.total_area)
    @test sys.lipid_weight == [1.0]
    @test sys.solvent_weight == [1.0]
    @test sys.flexibility == 0.25
    @test sys.lipid_head_tail_length[1] ≈ 11.574471473693848u"Å" rtol = 1e-6

    # invalid type
    @test_throws ArgumentError Membrane(
        type=:trilayer,
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    # mismatched lengths
    @test_throws ArgumentError Membrane(
        lipids=[lipid, lipid], lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    # neither / both of total_area, total_lipids given
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20, total_area=2000.0,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    # out-of-range / degenerate head-tail indices
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[9999], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[1], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    # flexibility out of (0, 0.5)
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20, flexibility=0.5,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )
    @test_throws ArgumentError Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20, flexibility=0.0,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    # density without units warns and assumes g/mL
    sys2 = @test_logs (:warn,) Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )
    @test sys2.solvent_density == 1.0u"g/mL"
end

@testitem "Membrane show method" begin
    using Packmol
    using Unitful
    using ShowMethodTesting

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    sys = Membrane(
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0u"g/mL",
    )
    @test parse_show(sys) ≈ """
    ==================================================================
    Membrane properties (bilayer):
    ==================================================================
        Lipids: palmitoil.pdb
        Lipid molar ratio (normalized): [1.0]
        Lipid head-to-tail lengths: 11.574471473693848 Å
        Solvents: water.pdb
        Solvent molar ratio (normalized): [1.0]
        Solvent density: 1.0 g mL^-1
        Solvent layer width: 15.0 Å
        Area per lipid: 60.0 Å^2
        Total lipids: 20
        Flexibility: 0.25
    ==================================================================
    """
end

@testitem "Membrane bilayer geometry" begin
    using Packmol
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    sys = Membrane(
        type=:bilayer,
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    tmp_input_file = tempname() * ".inp"
    rm(tmp_input_file, force=true)
    lipid_placements, solvent_placements, unitcell = write_packmol_input(
        sys; input=tmp_input_file, debug=true,
    )
    @test isfile(tmp_input_file)

    a, b, c, α, β, γ = Packmol._unitcell_abc_angles(unitcell)
    d = sys.lipid_head_tail_length[1]
    @test a ≈ b
    @test c ≈ ustrip(u"Å", 2 * 15.0u"Å" + 2d) # 2 solvent layers + 2 leaflets
    @test α ≈ β ≈ γ ≈ 90.0

    @test length(lipid_placements) == 2 # 1 lipid type x 2 leaflets
    @test length(solvent_placements) == 2 # 1 solvent x 2 slabs
    @test sum(lp.number for lp in lipid_placements) == 20
    @test all(sp.number > 0 for sp in solvent_placements)

    # leaflets are mirror images of each other about z = 0
    leaflet1, leaflet2 = lipid_placements
    @test leaflet1.zlo ≈ -leaflet2.zhi
    @test leaflet1.zhi ≈ -leaflet2.zlo
    @test leaflet1.head_low != leaflet2.head_low

    # The generated input file must be valid Packmol syntax for the native engine.
    psys = Packmol.read_packmol_input(tmp_input_file)
    @test psys.nmols == sum(lp.number for lp in lipid_placements) + sum(sp.number for sp in solvent_placements)
    rm(tmp_input_file, force=true)
end

@testitem "Membrane monolayer geometry" begin
    using Packmol
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    sys = Membrane(
        type=:monolayer,
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=10,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    tmp_input_file = tempname() * ".inp"
    rm(tmp_input_file, force=true)
    lipid_placements, solvent_placements, unitcell = write_packmol_input(
        sys; input=tmp_input_file, debug=true,
    )
    rm(tmp_input_file, force=true)

    a, b, c, = Packmol._unitcell_abc_angles(unitcell)
    d = sys.lipid_head_tail_length[1]
    @test c ≈ ustrip(u"Å", 15.0u"Å" + d) # 1 solvent layer + 1 leaflet

    @test length(lipid_placements) == 1
    @test length(solvent_placements) == 1
    leaflet = only(lipid_placements)
    slab = only(solvent_placements)
    # solvent is entirely above the leaflet, i.e. above the heads
    @test slab.zlo ≈ leaflet.zhi
    # the head sits at the high edge of the leaflet (facing the solvent)
    @test leaflet.head_low == false
end

@testitem "Membrane multi-lipid, multi-solvent, total_area" begin
    using Packmol
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    diatomic = joinpath(Packmol.src_dir, "..", "test", "structure_files", "diatomic.pdb")

    sys = Membrane(
        lipids=[lipid, lipid], lipid_head=[31, 30], lipid_tail=[1, 2], lipid_molar_ratio=[2.0, 1.0],
        area_per_lipid=60.0, total_area=2000.0,
        solvent=[water, diatomic], solvent_molar_ratio=[10.0, 1.0],
        solvent_layer_width=12.0, solvent_density=1.0,
    )
    @test sys.lipid_weight ≈ [2 / 3, 1 / 3]
    @test sys.solvent_weight ≈ [10 / 11, 1 / 11]

    tmp_input_file = tempname() * ".inp"
    rm(tmp_input_file, force=true)
    lipid_placements, solvent_placements, unitcell = write_packmol_input(
        sys; input=tmp_input_file, debug=true,
    )
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
    using PDBTools
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")

    sys = Membrane(
        type=:bilayer,
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=20,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    tmp_output_file = tempname() * ".pdb"
    packed_sys = packmol(sys; output=tmp_output_file, seed=42, iprint=0)
    @test packed_sys isa Packmol.PackmolSystem
    @test packed_sys.status == :packing_ready
    @test isfile(tmp_output_file)

    # Physical sanity: heads must end up near the outer (solvent-facing) edges
    # of the membrane and tails near the midplane, on both leaflets.
    atoms = read_pdb(tmp_output_file)
    lipid_natoms = length(read_pdb(lipid))
    n_lipids = 20
    d = ustrip(u"Å", sys.lipid_head_tail_length[1])
    flex = sys.flexibility * d
    for i in 1:n_lipids
        mol = atoms[(i-1)*lipid_natoms+1:i*lipid_natoms]
        head_z, tail_z = mol[31].z, mol[1].z
        @test abs(abs(head_z) - d) < flex + 3.0 # head near an outer edge (±d)
        @test abs(tail_z) < flex + 3.0          # tail near the midplane (z=0)
    end
    rm(tmp_output_file, force=true)
end

@testitem "Membrane packmol (monolayer)" begin
    using Packmol
    using PDBTools
    using Unitful

    lipid = joinpath(Packmol.src_dir, "..", "test", "structure_files", "palmitoil.pdb")
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")

    sys = Membrane(
        type=:monolayer,
        lipids=lipid, lipid_head=[31], lipid_tail=[1], lipid_molar_ratio=[1.0],
        area_per_lipid=60.0, total_lipids=10,
        solvent=water, solvent_layer_width=15.0, solvent_density=1.0,
    )

    tmp_output_file = tempname() * ".pdb"
    packed_sys = packmol(sys; output=tmp_output_file, seed=42, iprint=0)
    @test packed_sys isa Packmol.PackmolSystem
    @test packed_sys.status == :packing_ready
    @test isfile(tmp_output_file)
    rm(tmp_output_file, force=true)
end
