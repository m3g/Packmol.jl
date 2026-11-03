#
# Scaling of the rotations of the molecules in the packing loop.
#
# The atoms of a molecule are at xᵢ = cm + R(a)·rᵢ, with R(a) = Rx(β)·Ry(γ)·Rz(θ)
# (`eulermat`) and the reference coordinates rᵢ centered at the origin. For a change δa
# of the Euler angles, the atoms move by
#
#     Σᵢ |δxᵢ|² = δaᵀ·G·δa,   G = Jᵀ·I·J,
#
# where I = Σᵢ (|rᵢ|²·𝟙 - rᵢ·rᵢᵀ) is the inertia tensor of the reference coordinates (unit
# masses) and J maps the rates of the angles to the angular velocity in the frame of the
# molecule:
#
#     ω = J·(β̇, γ̇, θ̇),   J = [ cos γ cos θ   sin θ   0
#                               -cos γ sin θ   cos θ   0
#                                sin γ         0       1 ].
#
# The Euler angles are thus badly scaled variables: a change of an angle moves the atoms of
# a long molecule (a lipid) much more than those of a small one (a water molecule), and by
# an amount that depends on its current orientation. The solver has a single step length
# for all the variables (the spectral step of SPGBox), so it cannot correct for that.
#
# In each chunk of the packing loop, the angles of each molecule are replaced by
# a′ = Lᵀ·a, where G/N = L·Lᵀ (N atoms) at the angles of the start of the chunk: a unit
# change of a′ moves the atoms of the molecule by 1 Å RMS, as a unit change of the center
# of mass does. This is the Gauss-Newton approximation of the Hessian of the molecule,
# when all its atoms have the same curvature. The centers of mass, and thus their bounds,
# are unchanged. G is singular at γ = ±π/2 (gimbal lock), and for linear or one-atom
# molecules: it is regularized by a fraction of its trace, and replaced by the identity
# when it vanishes.
#

# Inertia tensor per atom of the (centered) reference coordinates
function _inertia_per_atom(ref::Vector{SVector{3,T}}) where {T}
    In = zero(SMatrix{3,3,T,9})
    for r in ref
        In += (r' * r) * one(SMatrix{3,3,T,9}) - r * r'
    end
    return In / length(ref)
end

# Lower Cholesky factor of a symmetric positive definite 3×3 matrix
function _cholesky_lower(A::SMatrix{3,3,T}) where {T}
    l11 = sqrt(A[1, 1])
    l21 = A[2, 1] / l11
    l31 = A[3, 1] / l11
    l22 = sqrt(A[2, 2] - l21^2)
    l32 = (A[3, 2] - l31 * l21) / l22
    l33 = sqrt(A[3, 3] - l31^2 - l32^2)
    z = zero(T)
    return SMatrix{3,3,T,9}(l11, l21, l31, z, l22, l32, z, z, l33)
end

# Factor L of the metric (per atom) of the Euler angles `angles`, for a molecule with
# inertia tensor per atom `Inp`: G/N = L·Lᵀ (regularized)
function _rotation_metric_factor(Inp::SMatrix{3,3,T}, angles::SVector{3,T}; regularization::T=T(1e-2)) where {T}
    _, γ, θ = angles
    c2, s2, c3, s3 = cos(γ), sin(γ), cos(θ), sin(θ)
    J = SMatrix{3,3,T,9}(c2 * c3, -c2 * s3, s2, s3, c3, zero(T), zero(T), zero(T), one(T))
    G = J' * Inp * J
    G = (G + G') / 2
    trG = G[1, 1] + G[2, 2] + G[3, 3]
    trG > eps(T) || return one(SMatrix{3,3,T,9})
    return _cholesky_lower(G + (regularization * trG / 3) * one(SMatrix{3,3,T,9}))
end

#
# Per free molecule: the inertia tensor per atom of its structure type, and the factors
# of the current chunk (Lᵀ, L⁻ᵀ, L⁻¹).
#
struct RotationScaling{T}
    inertia::Vector{SMatrix{3,3,T,9}}
    scaled::Vector{Bool}
    Lt::Vector{SMatrix{3,3,T,9}}
    Lt_inv::Vector{SMatrix{3,3,T,9}}
    L_inv::Vector{SMatrix{3,3,T,9}}
end

function RotationScaling(packmol_system::PackmolSystem{3,T}) where {T}
    inertia = SMatrix{3,3,T,9}[]
    scaled = Bool[]
    for st in packmol_system.structure_types
        st.fixed.fixed && continue
        Inp = _inertia_per_atom(st.reference_coordinates)
        scale = _scales_rotations(st.reference_coordinates, Inp)
        for _ in 1:st.number_of_molecules
            push!(inertia, Inp)
            push!(scaled, scale)
        end
    end
    n = length(inertia)
    id = one(SMatrix{3,3,T,9})
    return RotationScaling{T}(inertia, scaled, fill(id, n), fill(id, n), fill(id, n))
end

# Whether the rotations of a molecule are scaled: not for linear (or one-atom) molecules,
# whose metric is singular. A rotation about the axis of a linear molecule moves none of
# its atoms, so that a unit step of the scaled angles along it would be an arbitrarily
# large rotation, changing the orientation (and thus the metric) of the molecule far from
# the one the factors of the chunk were computed at.
function _scales_rotations(ref::Vector{SVector{3,T}}, Inp::SMatrix{3,3,T}) where {T}
    # The reference coordinates are aligned to their principal axes (see
    # `set_reference_coordinates!`), so the diagonal of the inertia tensor holds the
    # principal moments
    moments = (Inp[1, 1], Inp[2, 2], Inp[3, 3])
    return minimum(moments) >= T(0.05) * maximum(moments) && maximum(moments) > eps(T)
end

# Factors at the angles of `x_mol` (the free molecules, unscaled)
function _update_rotation_scaling!(rs::RotationScaling{T}, x_mol::AbstractVector{MoleculePosition{3,T}}) where {T}
    for k in eachindex(x_mol, rs.inertia)
        rs.scaled[k] || continue
        L = _rotation_metric_factor(rs.inertia[k], x_mol[k].angles)
        Lt = L'
        rs.Lt[k] = Lt
        rs.Lt_inv[k] = inv(Lt)
        rs.L_inv[k] = inv(L)
    end
    return rs
end

# a → a′ = Lᵀ·a, in place
function _scale_angles!(x_mol::AbstractVector{MoleculePosition{3,T}}, rs::RotationScaling{T}) where {T}
    for k in eachindex(x_mol)
        x_mol[k] = MoleculePosition(x_mol[k].cm, rs.Lt[k] * x_mol[k].angles)
    end
    return x_mol
end

# a′ → a = L⁻ᵀ·a′: `x_mol` from the scaled `y_mol`
function _unscale_angles!(
    x_mol::AbstractVector{MoleculePosition{3,T}}, y_mol::AbstractVector{MoleculePosition{3,T}}, rs::RotationScaling{T},
) where {T}
    for k in eachindex(x_mol, y_mol)
        x_mol[k] = MoleculePosition(y_mol[k].cm, rs.Lt_inv[k] * y_mol[k].angles)
    end
    return x_mol
end

# Gradient relative to the scaled angles: ∂f/∂a′ = L⁻¹·∂f/∂a, in place
function _scale_angle_gradient!(g_mol::AbstractVector{MoleculePosition{3,T}}, rs::RotationScaling{T}) where {T}
    for k in eachindex(g_mol)
        g_mol[k] = MoleculePosition(g_mol[k].cm, rs.L_inv[k] * g_mol[k].angles)
    end
    return g_mol
end
