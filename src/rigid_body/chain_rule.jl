#=

This functions apply the chain rule to compute the gradient of the minimum-distance
function relative to the rotations and translation of the rigid-body molecules, from 
the gradient computed from the displacement of the atoms in cartesian coordinates. 

This was implemented originally in the Fortran code of Packmol, by J. M. Martínez.

>> The function modifies the fg.g vector

=#
function chain_rule!(fg, packmol_system::PackmolSystem{D,T}) where {D,T}
    imol_offset = 0
    iat_offset = 0
    for structure_type in packmol_system.structure_types
        natoms_st = structure_type.natoms
        nmols_st = structure_type.number_of_molecules
        st_imol_offset = imol_offset
        st_iat_offset = iat_offset
        @sync for irange in chunks(1:nmols_st; n=Threads.nthreads())
            @spawn for i in irange
                imol = st_imol_offset + i
                ifmol = st_iat_offset + (i - 1) * natoms_st + 1
                ilmol = ifmol + natoms_st - 1
                gmol = partial_derivatives(
                    fg.g[imol],
                    packmol_system.molecule_positions[imol].angles,
                    structure_type.reference_coordinates,
                    @view(fg.gxcar[ifmol:ilmol]),
                )
                fg.g[imol] = gmol
            end
        end
        imol_offset += nmols_st
        iat_offset += nmols_st * natoms_st
    end
    return packmol_system
end

# Same as chain_rule! above, restricted to an explicit, possibly
# non-contiguous list of molecules. Used by adjust_constraints!'s active-set
# loop, alongside compute_positions_and_constraints_for_mols! — only
# molecules with a fresh fg.gxcar this pass need their gradient chain-ruled.
#
# partial_derivatives accumulates onto its `gmol` argument (see gcm/grot
# below), so chain_rule! normally relies on the caller having zeroed the
# *entire* fg.g array first (CellListMap.reset_output! in the full-system
# path). Here that would reintroduce an O(nmols) cost on every call, which is
# exactly what the active-set restriction is meant to avoid — so this passes
# a fresh zero seed per molecule directly instead of reading (possibly stale)
# fg.g[imol].
function chain_rule_for_mols!(
    fg, packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int}, mol_structure_type::Vector{Int}, mol_iat_first::Vector{Int},
) where {D,T}
    @sync for (_, mrange) in enumerate(chunks(1:length(mol_list); n=Threads.nthreads()))
        @spawn for k in mrange
            imol = mol_list[k]
            structure_type = packmol_system.structure_types[mol_structure_type[imol]]
            ifmol = mol_iat_first[imol]
            ilmol = ifmol + structure_type.natoms - 1
            gmol = partial_derivatives(
                zero(MoleculePosition{D,T}),
                packmol_system.molecule_positions[imol].angles,
                structure_type.reference_coordinates,
                @view(fg.gxcar[ifmol:ilmol]),
            )
            fg.g[imol] = gmol
        end
    end
    return packmol_system
end

#=

For a single molecule, computes the partial derivatives of the gradient of the
minimum distance between atoms relative to the rotations and translations of the
molecule, given the gradient of the minimum distance between atoms relative to the
displacement of the atoms in cartesian coordinates.

This function does not mutate any argument. It returns the gradient for the 
molecule as a MoleculePosition data structure, to be used in the chain_rule! function.

=#
function partial_derivatives(
    gmol::MoleculePosition{D,T}, 
    angles::SVector{D,T},
    reference_coordinates::Vector{SVector{D,T}},
    gxcar::AbstractVector{SVector{D,T}},
) where {D,T}
    # Compute rotation derivative matrices
    ∂v_matrices = rotation_derivative_matrices(angles)
    
    # Initialize gradients with current values
    gcm = gmol.cm
    grot = gmol.angles
    
    # Loop over all atoms in the molecule
    for i in eachindex(reference_coordinates, gxcar)
        x = reference_coordinates[i]
        gx = gxcar[i]
        # Gradient contribution from center of mass translation
        gcm += gx
        # Gradient contribution from rotations (chain rule)
        # dx_atom/d_angle_d = (∂R/∂angle_d) * x_ref
        # df/d_angle_d = sum_i dot(gx_car_i, (∂R/∂angle_d) * x_ref_i)
        grot += SVector(ntuple(d -> dot(gx, ∂v_matrices[d] * x), D))
    end
    return MoleculePosition{D,T}(gcm, grot)
end

#=

Compute the derivative matrices of the rotation matrix with respect to each angle.
Returns a tuple of D matrices, where D is the dimension (2 or 3).

For 2D: Returns (∂R/∂θ,) where θ is the single rotation angle
For 3D: Returns (∂R/∂β, ∂R/∂γ, ∂R/∂θ) for the three Euler angles

=#
function rotation_derivative_matrices(angles::SVector{3,T}) where {T}
    # 3D case: Euler angles (β, γ, θ)
    # Must match eulermat(β, γ, θ) from rigid_body.jl exactly.
    # eulermat uses: c1=cos(β), s1=sin(β), c2=cos(γ), s2=sin(γ), c3=cos(θ), s3=sin(θ)
    # Here: sb=s1, cb=c1, sg=s2, cg=c2, st=s3, ct=c3
    (sb, cb), (sg, cg), (st, ct) = sincos.(angles)
    #!format: off
    # ∂R/∂β: only rows 2,3 are nonzero (row 1 has no β dependence)
    ∂v∂β = @SMatrix[
             zero(T)            zero(T)   zero(T)
        -sb*st+cb*ct*sg   -sb*ct-cb*sg*st   -cb*cg
         cb*st+sb*ct*sg   -sb*sg*st+cb*ct   -sb*cg
    ]
    # ∂R/∂γ: all rows have γ dependence
    ∂v∂γ = @SMatrix[
        -sg*ct        sg*st      cg
         sb*cg*ct    -sb*cg*st   sb*sg
        -cb*cg*ct     cb*cg*st  -cb*sg
    ]
    # ∂R/∂θ: column 3 is zero (no θ dependence in column 3 of eulermat)
    ∂v∂θ = @SMatrix[
              -cg*st            -cg*ct    zero(T)
         cb*ct-sb*sg*st   -cb*st-sb*sg*ct zero(T)
         sb*ct+cb*sg*st    cb*sg*ct-sb*st  zero(T)
    ]
    #!format: on
    return (∂v∂β, ∂v∂γ, ∂v∂θ)
end

function rotation_derivative_matrices(angles::SVector{2,T}) where {T}
    # 2D case: Single rotation angle θ
    θ = angles[1]
    sθ, cθ = sincos(θ)
    #!format: off
    # Derivative of 2D rotation matrix with respect to θ
    # R(θ) = [cos(θ) -sin(θ)]    =>    dR/dθ = [-sin(θ) -cos(θ)]
    #        [sin(θ)  cos(θ)]                   [ cos(θ) -sin(θ)]
    ∂v∂θ = @SMatrix[
        -sθ  cθ
        -cθ -sθ
    ]
    #!format: on
    return (∂v∂θ,)
end

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
