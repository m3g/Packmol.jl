@testitem "constrain_rotation" setup=[TestFiles] begin
    using LinearAlgebra: norm, dot
    using StaticArrays: SVector
    using PDBToolsBase: position
    angle_z(v) = acosd(clamp(v[3] / norm(v), -1, 1))
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    # The vector from the center of atoms 1 and 2 to atom 42, within 20° of z
    st = structure_type(lipid; number=20,
        constraints=[InsideBox([-20, -20, -20], [20, 20, 20])],
        constrain_rotation=Dict(:z => (0.0, 20.0, [1, 2], [42])),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    packmol(sys)
    @test sys.status == :packing_ready
    x = [SVector{3,Float64}(position(at)) for at in get_atoms(sys)]
    angles = [angle_z(x[(i-1)*42+42] - (x[(i-1)*42+1] + x[(i-1)*42+2]) / 2) for i in 1:20]
    @test maximum(angles) < 21.0
    # Default vector: the longest axis of the molecule, here around 90° from x
    st = structure_type(lipid; number=20,
        constraints=[InsideBox([-20, -20, -20], [20, 20, 20])],
        constrain_rotation=Dict(:x => (90.0, 10.0)),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    packmol(sys)
    @test sys.status == :packing_ready
    rc = only(st.rotation_constraints)
    for mp in sys.molecule_positions
        v = Packmol.eulermat(mp.angles) * rc.reference_vector
        @test 79.0 < acosd(clamp(v[1] / norm(v), -1, 1)) < 101.0
    end
end

@testitem "rotation constraints: penalty, gradient and weights" setup=[TestFiles] begin
    using StaticArrays: SVector
    using LinearAlgebra: norm
    RC = Packmol.RotationConstraint{3,Float64}
    rc = RC(:z, 0.0, 10.0; from=[1], to=[2])
    @test rc.lower == 0.0 && rc.upper ≈ deg2rad(10)
    @test Packmol.rotation_penalty(rc, SVector(0.0, 0.1, 1.0)) == 0.0
    v = SVector(0.3, -0.5, 0.8)
    @test Packmol.rotation_penalty(rc, v) > 0
    # Gradient against finite differences
    h = 1e-6
    g = Packmol.rotation_gradient(rc, v)
    gnum = SVector(ntuple(i -> (Packmol.rotation_penalty(rc, v + h * SVector(ntuple(j -> j == i, 3))) -
        Packmol.rotation_penalty(rc, v - h * SVector(ntuple(j -> j == i, 3)))) / 2h, 3))
    # (at constant |v|, as for a rigid molecule: the component along v is not used)
    vhat = v / norm(v)
    @test g ≈ gnum - (gnum' * vhat) * vhat rtol = 1e-5
    # Lower bound (angle with x between 80° and 100°)
    rc = RC(:x, 90.0, 10.0)
    @test Packmol.rotation_penalty(rc, SVector(1.0, 0.0, 0.2)) > 0
    @test Packmol.rotation_penalty(rc, SVector(0.1, 0.0, 1.0)) == 0
    # Default weights: v = Σ wᵢ rᵢ is the z axis of the reference orientation, with the
    # length of the molecule along it
    ref = [SVector(cos(i), 0.5 * sin(i), 3.0 * i) for i in 1:10]
    Packmol.set_reference_coordinates!(ref)
    Packmol._set_rotation_constraint!(rc, ref)
    Lz = maximum(r -> r[3], ref) - minimum(r -> r[3], ref)
    @test rc.reference_vector ≈ SVector(0.0, 0.0, Lz)
    @test sum(rc.weights) ≈ 0 atol = 1e-12
    # Atom sets: the vector between their centers
    rc = RC(:z, 0.0, 10.0; from=[1, 2], to=[10])
    Packmol._set_rotation_constraint!(rc, ref)
    @test rc.reference_vector ≈ ref[10] - (ref[1] + ref[2]) / 2
    @test_throws ArgumentError RC(:w, 0.0, 10.0)
    @test_throws ArgumentError RC(:z, 0.0, 10.0; from=[1])
end

@testitem "rotation constraints: gradient of the objective function" setup=[TestFiles] begin
    using Packmol: InteratomicDistanceFG, compute_atom_positions!, initialize_molecules!, fg!
    using StaticArrays: SVector
    using CellListMap: ParticleSystem
    using Random
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    st = structure_type(lipid; number=4,
        constraints=[InsideBox([-10, -10, -10], [10, 10, 10])],
        constrain_rotation=Dict(:z => (0.0, 5.0, [1, 2], [40, 42]), :x => (90.0, 5.0)),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    initialize_molecules!(sys, Random.Xoshiro(1))
    ap = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    compute_atom_positions!(ap, sys.molecule_positions, sys)
    cl = ParticleSystem(xpositions=copy(ap), unitcell=[100.0, 100.0, 100.0], cutoff=2.2,
        output=InteratomicDistanceFG{3,Float64}(sys), output_name=:fg, parallel=false)
    free = collect(1:sys.nmols)
    x = copy(reinterpret(Float64, sys.molecule_positions))
    g = zeros(length(x))
    fg!(g, x, cl, sys, ap, free)
    @test sys.structure_types[1].rotation_constraints[1].reference_vector != zero(SVector{3,Float64})
    h = 1e-6
    gnum = map(eachindex(x)) do i
        xp = copy(x); xp[i] += h
        xm = copy(x); xm[i] -= h
        (fg!(zeros(length(x)), xp, cl, sys, ap, free) - fg!(zeros(length(x)), xm, cl, sys, ap, free)) / 2h
    end
    @test g ≈ gnum rtol = 1e-4
end
