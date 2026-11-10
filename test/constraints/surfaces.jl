@testsnippet SphereSurface begin
    using StaticArrays, LinearAlgebra
    # Points of a sphere of radius R (Fibonacci lattice), with the outward normals
    function sphere_points(R, n)
        pts = SVector{3,Float64}[]
        for i in 1:n
            θ = acos(1 - 2 * (i - 0.5) / n)
            φ = π * (1 + sqrt(5)) * i
            push!(pts, R * SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)))
        end
        return pts
    end
    fd_gradient(f, x; e = 1e-6) = SVector(ntuple(i -> (f(x + e * SVector(ntuple(j -> j == i ? 1.0 : 0.0, 3))) -
        f(x - e * SVector(ntuple(j -> j == i ? 1.0 : 0.0, 3)))) / 2e, 3))
end

@testitem "PointSurface distance" setup = [SphereSurface] begin
    pts = sphere_points(50.0, 3000)
    s = PointSurface(pts, pts)
    @test length(s.points) == 3000
    @test s.radius > 0
    @test occursin("3000 points with normals", sprint(show, MIME"text/plain"(), s))
    # The signed distance, near (smoothed: within the curvature bias) and far from the surface
    for (r, tol) in ((40.0, 0.3), (48.0, 0.3), (50.0, 0.3), (53.0, 0.3), (70.0, 0.5), (500.0, 1.0), (0.5, 1.0))
        for u in (SVector(1.0, 0.2, 0.1), SVector(-0.3, 0.8, -0.5), SVector(0.1, -0.2, -1.0))
            x = r * normalize(u)
            h, g = surface_distance(s, x)
            @test abs(h - (r - 50.0)) < tol
            @test dot(g, normalize(u)) > 0.9
        end
    end
    # Gradient of h, near the surface (where it is smooth)
    for x in (SVector(30.0, 25.0, 10.0), SVector(-20.0, 5.0, 44.0), SVector(3.0, -51.0, 4.0))
        h, g = surface_distance(s, x)
        @test isapprox(g, fd_gradient(y -> surface_distance(s, y)[1], x); atol = 1e-5)
    end
    # Without normals: the distance to the nearest point
    s0 = PointSurface(pts)
    for x in (SVector(10.0, 3.0, 2.0), SVector(60.0, -3.0, 2.0), SVector(400.0, 300.0, -100.0))
        d, g = surface_distance(s0, x)
        @test d ≈ minimum(norm(x - p) for p in pts)
        @test g ≈ normalize(x - pts[argmin([norm(x - p) for p in pts])])
    end
    @test_throws ArgumentError PointSurface(SVector{3,Float64}[])
    @test_throws ArgumentError PointSurface(pts, pts[1:10])
    @test_throws ArgumentError WithinSurface(s0, 2.0; offset = 1.0)
    @test_throws ArgumentError WithinSurface(s, -2.0)
end

@testitem "Surface constraints: penalty and gradient" setup = [SphereSurface] begin
    pts = sphere_points(50.0, 3000)
    s = PointSurface(pts, pts)
    within = WithinSurface(s, 3.0; offset = 10.0)    # 57 ≤ |x| ≤ 63
    notwithin = NotWithinSurface(s, 5.0)             # |x| ≤ 45 or |x| ≥ 55
    @test occursin("WithinSurface", sprint(show, MIME"text/plain"(), within))
    @test occursin("NotWithinSurface", sprint(show, MIME"text/plain"(), notwithin))
    u = normalize(SVector(0.3, -0.4, 0.8))
    for r in (60.0, 62.0)
        @test Packmol.constraint_penalty(within, r * u) == 0
        @test Packmol.constraint_gradient(within, r * u) == zero(u)
    end
    for r in (20.0, 45.0, 54.0, 66.0, 200.0)
        x = r * u
        @test Packmol.constraint_penalty(within, x) > 0
        @test isapprox(Packmol.constraint_gradient(within, x), fd_gradient(y -> Packmol.constraint_penalty(within, y), x); rtol = 1e-4)
    end
    for r in (10.0, 40.0, 58.0, 300.0)
        @test Packmol.constraint_penalty(notwithin, r * u) == 0
        @test Packmol.constraint_gradient(notwithin, r * u) == zero(u)
    end
    for r in (47.0, 49.0, 51.5, 54.0)
        x = r * u
        @test Packmol.constraint_penalty(notwithin, x) > 0
        @test isapprox(Packmol.constraint_gradient(notwithin, x), fd_gradient(y -> Packmol.constraint_penalty(notwithin, y), x); rtol = 1e-4)
    end
    # not within, beyond the smoothing radius of every point: satisfied without the far search
    @test Packmol._far_satisfies(notwithin)
    @test Packmol.constraint_penalty(notwithin, SVector(0.0, 0.0, 0.0)) == 0
end

@testitem "read_surface" setup = [SphereSurface] begin
    dir = mktempdir()
    pts = sphere_points(20.0, 500)
    # x y z nx ny nz (Å)
    file = joinpath(dir, "sphere.dat")
    open(file, "w") do io
        println(io, "# a sphere of radius 20")
        for p in pts
            println(io, join(p, " "), " ", join(normalize(p), " "))
        end
    end
    s = read_surface(file)
    @test length(s.points) == 500
    @test length(s.normals) == 500
    @test abs(surface_distance(s, SVector(25.0, 0.0, 0.0))[1] - 5.0) < 0.3
    # x y z
    file0 = joinpath(dir, "points.dat")
    open(io -> foreach(p -> println(io, join(p, " ")), pts), file0, "w")
    @test isempty(read_surface(file0).normals)
    @test read_surface(file0; radius = 7.0).radius == 7.0
    # A point file of TS2CG (nm)
    ts = joinpath(dir, "OuterBM.dat")
    open(ts, "w") do io
        println(io, "Box     100.000     100.000     100.000")
        println(io, "< Point NoPoints     $(length(pts))>")
        println(io, "< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2  >")
        println(io, "< Outer >")
        for (i, p) in enumerate(pts)
            n = normalize(p)
            println(io, i - 1, " 0 0.5 ", join(p ./ 10, " "), " ", join(n, " "), " 0 0 1 1 0 0 0.05 0.05")
        end
    end
    st = read_surface(ts)
    @test st.points ≈ s.points
    @test st.normals ≈ s.normals
    # Errors
    bad = joinpath(dir, "bad.dat")
    write(bad, "1.0 2.0\n")
    @test_throws ArgumentError read_surface(bad)
    write(bad, "1.0 2.0 3.0\n1.0 2.0 3.0 0.0 0.0 1.0\n")
    @test_throws ArgumentError read_surface(bad)
    write(bad, "1.0 2.0 x\n")
    @test_throws ArgumentError read_surface(bad)
    @test_throws ArgumentError read_surface(joinpath(dir, "missing.dat"))
end

@testitem "Surface constraints: input file and packing" setup = [TestFiles, SphereSurface] begin
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    cp(water, joinpath(dir, "water.pdb"))
    pts = sphere_points(15.0, 1500)
    mkdir(joinpath(dir, "surf"))
    open(joinpath(dir, "surf", "sphere.dat"), "w") do io
        for p in pts
            println(io, join(p, " "), " ", join(normalize(p), " "))
        end
    end
    # Water in a shell 4 Å thick, 5 Å outside the sphere; the oxygen atoms not within 4 Å
    # of the sphere itself (a redundant constraint, to test "not within")
    input = joinpath(dir, "shell.inp")
    write(input, """
        tolerance 2.0
        output shell.pdb
        structure water.pdb
          number 120
          within surface surf/sphere.dat 2.0 7.0
          atoms 1
            not within surface surf/sphere.dat 4.0
          end atoms
        end structure
        """)
    sys = redirect_stdout(devnull) do
        packmol(input)
    end
    @test sys.status == :packing_ready
    cs = sys.structure_types[1].constraints
    @test cs[1] isa Packmol.SurfaceDistance{Inside,Float64}
    @test cs[2] isa Packmol.SurfaceDistance{Outside,Float64}
    @test cs[1].surface === cs[2].surface   # the file is read once
    r = [norm(SVector(a.x, a.y, a.z)) for a in get_atoms(sys)]
    @test minimum(r) > 15.0 + 4.5
    @test maximum(r) < 15.0 + 9.5
    # Errors of the constraint line
    write(input, replace(read(input, String), "within surface surf/sphere.dat 2.0 7.0" => "within surface surf/sphere.dat"))
    @test_throws ArgumentError Packmol.read_packmol_input(input)
end
