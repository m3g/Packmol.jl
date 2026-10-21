@testitem "threaded outputs: sparse and dense reduction" begin
    using Packmol: InteratomicDistanceFG, read_packmol_input, _add_gxcar!, _add_fmol!
    import CellListMap
    using StaticArrays
    using Random

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
    ps = read_packmol_input(input_file)
    natoms, nmols = length(ps.atoms), ps.nmols
    main = InteratomicDistanceFG{3,Float64}(ps)
    rng = Random.MersenneTwister(1)
    # Few entries: sparse reduction; many entries: dense reduction
    for nentries in (5, 4 * natoms)
        copies = [CellListMap.copy_output(main) for _ in 1:4]
        expected_g = zeros(SVector{3,Float64}, natoms)
        expected_fmol = zeros(nmols)
        for _ in 1:nentries
            c = rand(rng, copies)
            i = rand(rng, 1:natoms)
            v = randn(rng, SVector{3,Float64})
            _add_gxcar!(c, i, v)
            expected_g[i] += v
            m = rand(rng, 1:nmols)
            w = rand(rng)
            _add_fmol!(c, m, w)
            expected_fmol[m] += w
        end
        CellListMap.reset_output!(main)
        CellListMap.reduce_output!(main, copies)
        @test main.gxcar ≈ expected_g
        @test main.fmol ≈ expected_fmol
        # The copies are left zeroed, ready for the next computation
        @test all(c -> all(iszero, c.gxcar) && all(iszero, c.fmol), copies)
        @test all(c -> isempty(c.touched_atoms) && isempty(c.touched_mols), copies)
    end
    # Resetting a copy zeroes the touched entries
    c = CellListMap.copy_output(main)
    _add_gxcar!(c, 3, SVector(1.0, 2.0, 3.0))
    _add_fmol!(c, 2, 1.0)
    CellListMap.reset_output!(c)
    @test all(iszero, c.gxcar) && all(iszero, c.fmol)
    @test isempty(c.touched_atoms) && isempty(c.touched_mols)
end

@testitem "constraint penalties per molecule" begin
    using Packmol: InteratomicDistanceFG, compute_atom_positions!, constraint_fg!, constraint_penalty
    using StaticArrays
    using Random
    # 101 waters (303 atoms): the chunks of atoms of the constraint computation
    # start and end in the middle of molecules
    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    water = structure_type(water_pdb; number=101, constraints=[InsideBox([0, 0, 0], [10, 10, 10])])
    ps = PackmolSystem([water]; output="", tolerance=2.0)
    Packmol._center_reference_coordinates!(ps)
    rng = Random.MersenneTwister(2)
    for imol in 1:ps.nmols
        ps.molecule_positions[imol] = Packmol.MoleculePosition(
            SVector{3}(30 .* rand(rng, 3) .- 10), SVector{3}(2π .* rand(rng, 3)))
    end
    x = Vector{SVector{3,Float64}}(undef, length(ps.atoms))
    compute_atom_positions!(x, ps.molecule_positions, ps)
    fg = InteratomicDistanceFG{3,Float64}(ps)
    constraint_fg!(fg, x, ps)
    c = ps.structure_types[1].constraints[1]
    fmol_ref = zeros(ps.nmols)
    for (iat, atom) in enumerate(ps.atoms)
        fmol_ref[atom.molecule_index] += constraint_penalty(c, x[iat])
    end
    @test sum(fmol_ref) > 0
    @test fg.fmol ≈ fmol_ref
    @test fg.f ≈ sum(fmol_ref)
end

@testitem "parallel and serial fg! agree" begin
    using Packmol: InteratomicDistanceFG, compute_bounding_box, compute_atom_positions!,
        initialize_molecules!, fg!, read_packmol_input
    using StaticArrays
    using CellListMap: ParticleSystem
    using Random

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
    ps = read_packmol_input(input_file)
    initialize_molecules!(ps, Random.MersenneTwister(42))
    atom_positions = Vector{SVector{3,Float64}}(undef, length(ps.atoms))
    compute_atom_positions!(atom_positions, ps.molecule_positions, ps)
    lo, hi = compute_bounding_box(atom_positions)
    free_mol_indices = collect(1:ps.nmols)
    x0 = copy(reinterpret(Float64, ps.molecule_positions))
    # Serial, and parallel with 4 batches (the threaded copies of the output are
    # used regardless of the number of threads)
    results = map(((false, (0, 0)), (true, (1, 4)))) do (parallel, nbatches)
        cl_system = ParticleSystem(
            xpositions=copy(atom_positions),
            unitcell=3 * (hi - lo),
            cutoff=2.2,
            output=InteratomicDistanceFG{3,Float64}(ps),
            output_name=:fg,
            parallel=parallel,
            nbatches=nbatches,
        )
        x = copy(x0)
        g = zeros(length(x))
        # two evaluations, such that the second one starts from the copies left
        # by the first one
        fg!(g, x, cl_system, ps, atom_positions, free_mol_indices)
        x .*= 1.01
        f = fg!(g, x, cl_system, ps, atom_positions, free_mol_indices)
        (f, g, copy(cl_system.fg.fmol))
    end
    (f1, g1, fmol1), (f2, g2, fmol2) = results
    @test f1 > 0
    @test f1 ≈ f2
    @test g1 ≈ g2
    @test fmol1 ≈ fmol2
end
