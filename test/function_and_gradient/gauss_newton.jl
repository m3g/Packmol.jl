@testitem "Gauss-Newton products" setup=[TestFiles] begin
    using Packmol: InteratomicDistanceFG, compute_bounding_box, compute_atom_positions!,
        initialize_molecules!, fg!, read_packmol_input, gauss_newton_product!, GaussNewtonBuffers,
        MoleculePosition, _build_mol_structure_type, _build_mol_iat_first
    using StaticArrays
    using LinearAlgebra: dot
    using CellListMap: ParticleSystem
    using Random

    input_file = TestFiles.test_file("input_files/water_box_small.inp")
    ps = read_packmol_input(input_file)
    initialize_molecules!(ps, Random.MersenneTwister(42))
    natoms = length(ps.atoms)
    atom_positions = Vector{SVector{3,Float64}}(undef, natoms)
    compute_atom_positions!(atom_positions, ps.molecule_positions, ps)
    lo, hi = compute_bounding_box(atom_positions)
    free_mol_indices = collect(1:ps.nmols)
    mol_structure_type = _build_mol_structure_type(ps)
    mol_iat_first = _build_mol_iat_first(ps)
    x0 = copy(reinterpret(Float64, ps.molecule_positions))
    n = length(x0)
    rng = Random.MersenneTwister(7)
    v = randn(rng, n)
    u = randn(rng, n)

    # Displacements of the atoms for a change v of the variables, by finite differences
    function atom_displacements(x, v; h=1e-6)
        pos(x) = (ap = similar(atom_positions); compute_atom_positions!(ap, reinterpret(MoleculePosition{3,Float64}, x), ps); ap)
        (pos(x .+ h .* v) .- pos(x .- h .* v)) ./ (2h)
    end

    # Serial, and parallel with 4 batches (the recorded terms are merged from the copies)
    results = map(((false, (0, 0)), (true, (1, 4)))) do (parallel, nbatches)
        output = InteratomicDistanceFG{3,Float64}(ps)
        output.record = true
        cl_system = ParticleSystem(
            xpositions=copy(atom_positions),
            unitcell=3 * (hi - lo),
            cutoff=2.2,
            output=output,
            output_name=:fg,
            parallel=parallel,
            nbatches=nbatches,
        )
        x = copy(x0)
        g = zeros(n)
        f = fg!(g, x, cl_system, ps, atom_positions, free_mol_indices)
        buffers = GaussNewtonBuffers{3,Float64}(natoms)
        Hv = zeros(n)
        Hu = zeros(n)
        gauss_newton_product!(Hv, x, v, cl_system.fg, ps, free_mol_indices, mol_structure_type, mol_iat_first, buffers)
        gauss_newton_product!(Hu, x, u, cl_system.fg, ps, free_mol_indices, mol_structure_type, mol_iat_first, buffers)
        # The quadratic form, from the recorded terms and the displacements of the atoms
        dx = atom_displacements(x, v)
        q = sum(2 * dot(p.u, dx[p.j] - dx[p.i])^2 for p in cl_system.fg.gn_pairs; init=0.0) +
            sum(dot(c.g, dx[c.iat])^2 / (2 * c.p) for c in cl_system.fg.gn_constraints; init=0.0)
        (; f, Hv, Hu, q, npairs=length(cl_system.fg.gn_pairs))
    end
    r1, r2 = results
    @test r1.npairs > 0
    @test r1.npairs == r2.npairs
    @test r1.Hv ≈ r2.Hv
    # vᵀHv is the quadratic form of the recorded terms, H is symmetric
    @test dot(v, r1.Hv) ≈ r1.q rtol = 1e-6
    @test dot(v, r1.Hv) > 0
    @test dot(u, r1.Hv) ≈ dot(v, r1.Hu)
end
