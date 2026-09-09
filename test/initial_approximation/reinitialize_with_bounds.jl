@testitem "reinitialize_with_bounds! doesn't inherit Step 2's edge-pileup" begin
    # Regression test for a real bug: a structure type confined to a thin PBC
    # slab by two one-sided `plane` constraints, where one of the two planes
    # coincides with the periodic wrap seam (so it's essentially never
    # violated from a uniform starting draw), used to end up with almost all
    # of its molecules piled up hugging the *other* plane instead of spread
    # uniformly across the slab. Root cause: `reinitialize_with_bounds!`
    # skipped re-randomizing any molecule whose Step-2 (constraint-only fit)
    # position already satisfied its own constraint precision — which
    # Step 2's own nearest-point gradient descent guarantees for exactly the
    # molecules that got pushed to that single edge, so the pileup survived
    # untouched instead of being redrawn uniformly. Reported in practice as a
    # wildly non-uniform solvent density profile in a `Membrane` bilayer
    # (dense near the membrane, empty near the box's outer PBC face) that
    # persisted even after `check=true`-only (no optimization loop at all).
    using Packmol
    using StaticArrays
    using LinearAlgebra

    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")

    Lx, Ly, Lz = 20.0, 20.0, 60.0
    w = 15.0 # slab width
    zlo = -Lz / 2       # coincides with the periodic wrap seam: essentially never violated
    zhi = -Lz / 2 + w   # the only edge Step 2's gradient descent ever actually pushes molecules toward

    st = structure_type(water; number=600, constraints=[AbovePlane([0, 0, 1], zlo), BelowPlane([0, 0, 1], zhi)])
    unitcell = Matrix{Float64}(Diagonal([Lx, Ly, Lz]))
    sys = PackmolSystem([st]; output="", unitcell, unitcell_center=zero(SVector{3,Float64}), check=true, seed=42)

    Packmol.packmol(sys)
    @test sys.status == :not_packed # check mode: initial approximation only, no optimization loop

    zs = [mp.cm[3] for mp in sys.molecule_positions]
    @test all(z -> zlo - 1e-6 <= z <= zhi + 1e-6, zs) # all within the slab (up to the confinement's own tolerance)

    # Split the slab into a half nearest the "always satisfied" edge (zlo, the
    # wrap seam) and a half nearest the "actually enforced" edge (zhi): a
    # uniform distribution puts about half the molecules in each. Before the
    # fix, this reproduced a ~9x skew toward the zhi half; a generous 2x
    # bound here still clearly catches a regression without being flaky.
    zmid = (zlo + zhi) / 2
    n_near_zlo = count(z -> z < zmid, zs)
    n_near_zhi = count(z -> z >= zmid, zs)
    @test n_near_zlo > 0
    @test n_near_zhi > 0
    @test max(n_near_zlo, n_near_zhi) / min(n_near_zlo, n_near_zhi) < 2.0
end
