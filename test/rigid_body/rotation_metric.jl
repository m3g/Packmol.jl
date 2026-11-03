@testitem "rotation metric" begin
    using Packmol: eulermat, _inertia_per_atom, _rotation_metric_factor, _cholesky_lower
    using StaticArrays
    using LinearAlgebra: norm, I
    using Random: Xoshiro

    rng = Xoshiro(1)
    # An elongated molecule, centered
    ref = [SVector{3}(randn(rng, 3) .* (1.0, 2.0, 8.0)) for _ in 1:40]
    c = sum(ref) / length(ref)
    ref = [r - c for r in ref]
    Inp = _inertia_per_atom(ref)

    # Cholesky factor
    A = SMatrix{3,3}(4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0)
    L = _cholesky_lower(A)
    @test L * L' ≈ A
    @test L[1, 2] == L[1, 3] == L[2, 3] == 0

    # The RMS displacement of the atoms for a small change δa of the angles is
    # |Lᵀ·δa| (up to the regularization, which is not used here)
    for _ in 1:10
        a = SVector{3}(2π .* rand(rng, 3) .- π)
        δ = SVector{3}(1e-6 .* randn(rng, 3))
        R0, R1 = eulermat(a), eulermat(a + δ)
        rms2 = sum(sum(abs2, (R1 - R0) * r) for r in ref) / length(ref)
        L = _rotation_metric_factor(Inp, a; regularization=0.0)
        @test isapprox(rms2, sum(abs2, L' * δ); rtol=1e-4)
    end

    # Near the gimbal lock (γ = π/2) and for a one-atom molecule, the factor is invertible
    L = _rotation_metric_factor(Inp, SVector(0.3, π / 2, 0.1))
    @test all(isfinite, inv(L))
    @test _rotation_metric_factor(_inertia_per_atom([SVector(0.0, 0.0, 0.0)]), SVector(0.1, 0.2, 0.3)) ≈ I
end

@testitem "rotation scaling of the variables" begin
    using Packmol: MoleculePosition, RotationScaling, _update_rotation_scaling!, _scale_angles!,
        _unscale_angles!, _scale_angle_gradient!, _inertia_per_atom, _rotation_metric_factor
    using StaticArrays
    using Random: Xoshiro

    rng = Xoshiro(2)
    ref = [SVector{3}(randn(rng, 3) .* (1.0, 1.5, 6.0)) for _ in 1:20]
    c = sum(ref) / length(ref)
    ref = [r - c for r in ref]
    Inp = _inertia_per_atom(ref)
    n = 5
    id = one(SMatrix{3,3,Float64,9})
    rs = RotationScaling{Float64}(fill(Inp, n), fill(id, n), fill(id, n), fill(id, n))
    x = [MoleculePosition(SVector{3}(randn(rng, 3)), SVector{3}(randn(rng, 3))) for _ in 1:n]
    _update_rotation_scaling!(rs, x)

    # Round trip, centers of mass unchanged
    y = copy(x)
    _scale_angles!(y, rs)
    @test all(y[k].cm == x[k].cm for k in 1:n)
    x2 = similar(x)
    _unscale_angles!(x2, y, rs)
    @test all(x2[k].cm ≈ x[k].cm && x2[k].angles ≈ x[k].angles for k in 1:n)

    # The gradient relative to the scaled angles is the derivative of a function of the
    # scaled angles: f(a) = Σ |R(a)·r|² weighted, compared with finite differences in a′
    w = SVector(0.3, -1.2, 0.7)
    f(xs) = sum(sum(w' * (Packmol.eulermat(m.angles) * r) for r in ref) for m in xs)
    function grad(xs)
        h = 1e-6
        [MoleculePosition(zero(SVector{3,Float64}), SVector{3}(ntuple(d -> begin
            e = SVector{3}(ntuple(i -> i == d ? h : 0.0, 3))
            xp = copy(xs); xp[k] = MoleculePosition(xs[k].cm, xs[k].angles + e)
            xm = copy(xs); xm[k] = MoleculePosition(xs[k].cm, xs[k].angles - e)
            (f(xp) - f(xm)) / 2h
        end, 3))) for k in eachindex(xs)]
    end
    g = grad(x)
    _scale_angle_gradient!(g, rs)
    fy(ys) = (xs = similar(ys); _unscale_angles!(xs, ys, rs); f(xs))
    h = 1e-6
    for k in 1:n, d in 1:3
        e = SVector{3}(ntuple(i -> i == d ? h : 0.0, 3))
        yp = copy(y); yp[k] = MoleculePosition(y[k].cm, y[k].angles + e)
        ym = copy(y); ym[k] = MoleculePosition(y[k].cm, y[k].angles - e)
        @test isapprox((fy(yp) - fy(ym)) / 2h, g[k].angles[d]; rtol=1e-4, atol=1e-7)
    end
end
