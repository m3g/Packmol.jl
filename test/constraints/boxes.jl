@testitem "Box constructors" begin
    # Positional: the parameters of the input file (lower corner and side; two corners)
    @test InsideCube([0,0,0],1.) == Cube{Inside,Float64}([0.5,0.5,0.5],1.,5.0)
    @test InsideCube([0,0,0],1.,2.0) == Cube{Inside,Float64}([0.5,0.5,0.5],1.,2.0)
    @test OutsideCube([0,0,0],1.) == Cube{Outside,Float64}([0.5,0.5,0.5],1.,5.0)
    @test InsideBox([0,0,0],[1.,2.,4.]) == Box{Inside,Float64}([0.5,1.,2.],[1.,2.,4.],5.0)
    @test OutsideBox([0,0,0],[1.,2.,4.]) == Box{Outside,Float64}([0.5,1.,2.],[1.,2.,4.],5.0)
    @test InsideBox([0,0,0],[1.,2.,4.]; weight=2.0) == Box{Inside,Float64}([0.5,1.,2.],[1.,2.,4.],2.0)
    @test_throws ArgumentError InsideBox([0,0,0],[1.,-1.,1.])
    @test_throws ArgumentError OutsideBox([0,0,0],[1.,-1.,1.])
    # Keywords: center and sides
    @test InsideCube(center=[0,0,0],side=1.) == Cube{Inside,Float64}([0.,0.,0.],1.,5.0)
    @test InsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Inside,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideCube(center=[0,0,0],side=1.) == Cube{Outside,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Outside,Float64}([0.,0.,0.],1.,2.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Inside,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Inside,Float64}([0.,0.,0.],[1.,1.,1.],2.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Outside,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Outside,Float64}([0.,0.,0.],[1.,1.,1.],2.0)
    # The input file lines give the same constraints as the positional constructors
    parse(line) = Packmol.parse_constraint(Val(Symbol(join(split(line)[1:2], " "))), split(line)[3:end], Float64, "")
    @test parse("inside box 0. 1. 2. 10. 11. 12.") == InsideBox([0., 1., 2.], [10., 11., 12.])
    @test parse("outside box 0. 1. 2. 10. 11. 12.") == OutsideBox([0., 1., 2.], [10., 11., 12.])
    # As in Fortran Packmol, `inside cube xmin ymin zmin side`
    @test parse("inside cube 0. 1. 2. 10.") == InsideCube([0., 1., 2.], 10.)
    @test parse("outside cube 0. 1. 2. 10.") == OutsideCube([0., 1., 2.], 10.)
    c = InsideCube([0., 0., 0.], 10.)
    @test Packmol.constraint_penalty(c, [9.5, 9.5, 9.5]) == 0
    @test Packmol.constraint_penalty(c, [-0.5, 5., 5.]) > 0
end

@testitem "Constraint gradients" begin
    using ForwardDiff
    using StaticArrays
    x = SVector{3,Float64}(1.5,1.0,0.)
    c = InsideCube([0.2, 0., 0.1], 0.5)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideCube([0.2, 0., 0.1], 2.)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = InsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
end

@testitem "Outside box and cube" begin
    using ForwardDiff
    using StaticArrays
    # Outside if outside in any dimension: a slab |z| < 24 excluded from a large box
    c = OutsideBox([-100., -100., -24.], [100., 100., 24.])
    @test Packmol.constraint_penalty(c, SVector(5.0, 3.0, 40.0)) == 0
    @test Packmol.constraint_penalty(c, SVector(5.0, 3.0, -30.0)) == 0
    @test Packmol.constraint_gradient(c, SVector(5.0, 3.0, 40.0)) == zeros(3)
    # Inside: the squared distance to the nearest face, pushed out through it
    @test Packmol.constraint_penalty(c, SVector(5.0, 3.0, 20.0)) ≈ 5 * 4^2
    @test Packmol.constraint_gradient(c, SVector(5.0, 3.0, 20.0)) ≈ [0, 0, -2 * 5 * 4]
    @test Packmol.constraint_gradient(c, SVector(5.0, 3.0, -20.0)) ≈ [0, 0, 2 * 5 * 4]
    cube = OutsideCube([0., 0., 0.], 10.)
    @test Packmol.constraint_penalty(cube, SVector(1.0, 5.0, 5.0)) ≈ 5 * 1^2
    @test Packmol.constraint_penalty(cube, SVector(11.0, 5.0, 5.0)) == 0
    for c in (cube, OutsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0]), OutsideBox([-1., -2., -3.], [1., 2., 3.]))
        for x in (SVector(0.25, 0.2, 0.5), SVector(-0.5, 1.1, 0.0), SVector(3.0, 0.1, 0.3), SVector(2.0, 6.0, 4.5))
            @test ForwardDiff.gradient(y -> Packmol.constraint_penalty(c, y), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
end
