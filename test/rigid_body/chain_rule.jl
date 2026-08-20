@testitem "rotation derivative matrices" begin
    using Packmol: rotation_derivative_matrices, eulermat
    using StaticArrays
    using LinearAlgebra: norm

    for angles in [
        SVector(0.3, 0.7, 1.1),
        SVector(0.0, 0.0, 0.0),
        SVector(π / 3, π / 4, π / 6),
        SVector(1.5, -0.5, 2.0),
    ]
        β, γ, θ = angles
        ∂β, ∂γ, ∂θ = rotation_derivative_matrices(angles)

        h = 1e-7
        num_∂β = (eulermat(β + h, γ, θ) - eulermat(β - h, γ, θ)) / (2h)
        num_∂γ = (eulermat(β, γ + h, θ) - eulermat(β, γ - h, θ)) / (2h)
        num_∂θ = (eulermat(β, γ, θ + h) - eulermat(β, γ, θ - h)) / (2h)

        @test norm(∂β - num_∂β) < 1e-6
        @test norm(∂γ - num_∂γ) < 1e-6
        @test norm(∂θ - num_∂θ) < 1e-6
    end
end

@testitem "gradient chain rule 3D" begin
    using Packmol
    using Packmol: InteratomicDistanceFG, compute_bounding_box, compute_atom_positions!,
        initialize_molecules!, fg!, read_packmol_input, MoleculePosition, adjust_constraints!
    using StaticArrays
    using LinearAlgebra: norm
    import CellListMap
    using CellListMap: ParticleSystem
    using Random

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
    ps = read_packmol_input(input_file; D=3, T=Float64)

    rng = Random.MersenneTwister(42)
    initialize_molecules!(ps, rng)

    D = 3; T = Float64
    natoms = length(ps.atoms)

    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, ps.molecule_positions, ps)

    lo, hi = compute_bounding_box(atom_positions)
    box_size = hi - lo
    unitcell = T(3) * box_size
    tol = ps.tolerance
    packing_tol = tol + tol / 10
    volume = prod(unitcell)
    ncells = min(volume / packing_tol^D, natoms)
    cutoff = max(packing_tol, (volume / ncells)^(one(T) / D))

    fg_output = InteratomicDistanceFG{D,T}(ps)
    cl_system = ParticleSystem(
        xpositions=atom_positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=fg_output,
        output_name=:fg,
        parallel=false,
    )

    free_mol_indices = collect(1:ps.nmols)

    x = copy(reinterpret(T, ps.molecule_positions))
    g = zeros(length(x))
    f = fg!(g, x, cl_system, ps, atom_positions, free_mol_indices)
    g_analytical = copy(g)

    # Numerical gradient via central finite differences
    h = 1e-5
    g_num = zeros(length(x))
    for i in eachindex(x)
        x_save = x[i]
        x[i] = x_save + h
        g_tmp = zeros(length(x))
        f_plus = fg!(g_tmp, x, cl_system, ps, atom_positions, free_mol_indices)
        x[i] = x_save - h
        f_minus = fg!(g_tmp, x, cl_system, ps, atom_positions, free_mol_indices)
        x[i] = x_save
        g_num[i] = (f_plus - f_minus) / (2h)
    end

    @test maximum(abs.(g_analytical - g_num)) < 1e-3
    @test norm(g_analytical) > 0  # nonzero gradients
    @test norm(g_analytical - g_num) / norm(g_num) < 1e-4
end

@testitem "gradient chain rule 2D" begin

end
