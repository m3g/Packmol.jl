#
# Products with the Gauss-Newton approximation of the Hessian of the objective function.
#
# The objective function is a sum of squares, f = Σₚ rₚ², of the residuals of the
# overlapping pairs of atoms, rₚ = dₚ - tolₚ, and of the violated constraints, whose
# penalties are taken as the squares of single residuals, rₖ = √pₖ. Its Hessian is
# 2·Σ (∇rₚ ∇rₚᵀ + rₚ ∇²rₚ). The Gauss-Newton approximation keeps the first term only,
# which is positive semi-definite, and exact where the residuals vanish (near a solution):
#
#     H ≈ 2·Kᵀ·(Σₚ aₚ aₚᵀ)·K,
#
# with aₚ the gradients of the residuals relative to the positions of the atoms and
# K = ∂x/∂q the map of the variables of the molecules (centers of mass and Euler angles)
# to the positions of their atoms. For a pair, aₚ is u (the unit vector from atom i to
# atom j) at atom j and -u at atom i; for a constraint, aₖ = ∇pₖ/(2√pₖ) at its atom, so
# that 2·aₖaₖᵀ = ∇pₖ∇pₖᵀ/(2pₖ).
#
# The terms (the pairs and the violated constraints, with their vectors) are the ones
# recorded by the last evaluation of `fg!` with `record` set (see GaussNewtonPair in
# interatomic_distance_fg.jl): the product is a loop over them, without a search of
# neighbors. Rotation constraints are not included.
#

# Buffers of the product: displacements of the atoms (δx = K·v) and the products at the
# atoms (w = Σ aaᵀ·δx), for all atoms (zero for the fixed ones)
struct GaussNewtonBuffers{D,T}
    dx::Vector{SVector{D,T}}
    w::Vector{SVector{D,T}}
end
GaussNewtonBuffers{D,T}(natoms::Int) where {D,T} =
    GaussNewtonBuffers{D,T}(fill(zero(SVector{D,T}), natoms), fill(zero(SVector{D,T}), natoms))

"""
    gauss_newton_product!(Hv, x, v, fg_output, packmol_system, free_mol_indices, mol_structure_type, mol_iat_first, buffers)

`Hv = H·v`, with H the Gauss-Newton approximation of the Hessian of the objective function
at the variables `x` (of the free molecules, with the Euler angles), from the terms recorded
in `fg_output` by the last evaluation of `fg!`.
"""
function gauss_newton_product!(
    Hv::AbstractVector{T}, x::AbstractVector{T}, v::AbstractVector{T},
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    free_mol_indices::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
    buffers::GaussNewtonBuffers{D,T},
) where {D,T}
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    v_mol = reinterpret(MoleculePosition{D,T}, v)
    Hv_mol = reinterpret(MoleculePosition{D,T}, Hv)
    (; dx, w) = buffers
    @_sync for irange in index_chunks(w; n=_nbatches())
        @spawn for i in irange
            w[i] = zero(SVector{D,T})
        end
    end
    # Displacements of the atoms of the free molecules: δxᵢ = δcm + (Σ_d ∂R/∂a_d δa_d)·rᵢ
    @_sync for krange in index_chunks(1:length(free_mol_indices); n=_nbatches())
        @spawn for k in krange
            imol = free_mol_indices[k]
            st = packmol_system.structure_types[mol_structure_type[imol]]
            ref = st.reference_coordinates
            ∂R = rotation_derivative_matrices(x_mol[k].angles)
            δa = v_mol[k].angles
            dR = sum(∂R[d] * δa[d] for d in 1:D)
            δcm = v_mol[k].cm
            iat_first = mol_iat_first[imol]
            for j in eachindex(ref)
                dx[iat_first+j-1] = δcm + dR * ref[j]
            end
        end
    end
    # w = Σ aaᵀ·δx (times 2), over the recorded pairs and constraints
    for pair in fg_output.gn_pairs
        u = pair.u
        s = 2 * dot(u, dx[pair.j] - dx[pair.i])
        w[pair.j] += s * u
        w[pair.i] -= s * u
    end
    for c in fg_output.gn_constraints
        w[c.iat] += (dot(c.g, dx[c.iat]) / (2 * c.p)) * c.g
    end
    # Back to the variables of the molecules: Kᵀ·w (the chain rule of the gradient)
    @_sync for krange in index_chunks(1:length(free_mol_indices); n=_nbatches())
        @spawn for k in krange
            imol = free_mol_indices[k]
            st = packmol_system.structure_types[mol_structure_type[imol]]
            iat_first = mol_iat_first[imol]
            Hv_mol[k] = partial_derivatives(
                zero(MoleculePosition{D,T}), x_mol[k].angles, st.reference_coordinates,
                @view(w[iat_first:iat_first+st.natoms-1]),
            )
        end
    end
    return Hv
end

# (Fixed atoms are never displaced: their entries of `dx` stay zero.)
