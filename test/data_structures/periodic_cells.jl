@testitem "Wigner-Seitz cells" begin
    using LinearAlgebra: norm, I
    # Truncated octahedron: 8 hexagonal faces at d/2 and 6 square faces at d/√3
    uc = octahedral_unitcell(Float64, 10.0)
    @test Packmol._is_wigner_seitz_cell(uc)
    faces = Packmol._wigner_seitz_face_vectors(uc)
    @test length(faces) == 14
    @test count(u -> norm(u) / 2 ≈ 5.0, faces) == 8
    @test count(u -> norm(u) / 2 ≈ 10 / sqrt(3), faces) == 6
    # Rhombic dodecahedron: 12 faces at d/2
    uc = Packmol.dodecahedral_unitcell(Float64, 10.0)
    @test Packmol._is_wigner_seitz_cell(uc)
    faces = Packmol._wigner_seitz_face_vectors(uc)
    @test length(faces) == 12
    @test all(u -> norm(u) / 2 ≈ 5.0, faces)
    @test Packmol._wigner_seitz_half_extent(uc) ≈ [5.0, 5.0, 5 * sqrt(2)]
    # Other cells are not
    @test !Packmol._is_wigner_seitz_cell(Matrix(10.0I, 3, 3))
    @test !Packmol._is_wigner_seitz_cell([10.0 5.0 0; 0 10 0; 0 0 10])
end

@testitem "constraints in the truncated octahedron" setup=[TestFiles] begin
    using LinearAlgebra: norm
    # Water outside a sphere larger than the sphere inscribed in the parallelepiped of
    # the cell vectors (0.408 d), but smaller than the one inscribed in the truncated
    # octahedron (d/2): the constraint holds in the truncated octahedron, centered at the
    # origin, and no water is placed inside the sphere through a periodic image
    d, r = 30.0, 13.5
    water = structure_type(TestFiles.test_file("structure_files/water.pdb");
        number = 300,
        constraints = [OutsideSphere([0.0, 0.0, 0.0], r)],
    )
    system = PackmolSystem([water];
        output = tempname() * ".pdb",
        tolerance = 2.0,
        unitcell = octahedral_unitcell(Float64, d),
    )
    @test system.wigner_seitz_cell
    redirect_stdout(devnull) do
        packmol(system)
    end
    @test system.status == :packing_ready
    triclinic_to_octahedral(system)
    atoms = get_atoms(system)
    @test minimum(norm([a.x, a.y, a.z]) for a in atoms) > r - 0.1
    # The whole truncated octahedron is filled: along z it extends to 0.611 d, beyond
    # the parallelepiped (0.408 d), and that region is only reached through periodic
    # images of the positions drawn in the parallelepiped
    @test maximum(abs(a.z) for a in atoms) > 0.55 * d
end
