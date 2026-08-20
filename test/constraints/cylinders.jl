@testitem "Cylinder constructors" begin
    @test InsideCylinder([0,0,0],[0,0,1],1.,2.) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test InsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test InsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.,weight=2.0) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,2.0)
    @test OutsideCylinder([0,0,0],[0,0,1],1.,2.) == Cylinder{Outside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test OutsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.) == Cylinder{Outside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
end

@testitem "Cylinder gradients" begin
    using ForwardDiff
    using StaticArrays
    # Axis-aligned cylinder along z, from z=0 to z=2, radius 1
    c = InsideCylinder(center=[0.0, 0.0, 0.0], axis=[0.0, 0.0, 1.0], radius=1.0, length=2.0)
    for x in (
        SVector(0.3, 0.2, 1.0),   # inside
        SVector(0.3, 0.2, -0.5),  # behind near cap
        SVector(0.3, 0.2, 2.7),   # beyond far cap
        SVector(1.5, 0.2, 1.0),   # outside radius
        SVector(1.5, 0.2, -0.5),  # behind cap and outside radius
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    end

    c2 = OutsideCylinder(center=[0.0, 0.0, 0.0], axis=[0.0, 0.0, 1.0], radius=1.0, length=2.0)
    for x in (
        SVector(0.3, 0.2, 1.0),   # inside the forbidden region
        SVector(0.3, 0.2, -0.5),  # behind near cap (already outside)
        SVector(0.3, 0.2, 2.7),   # beyond far cap (already outside)
        SVector(1.5, 0.2, 1.0),   # outside radius (already outside)
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c2, x), x) ≈ Packmol.constraint_gradient(c2, x)
    end

    # A tilted, off-center, unnormalized-axis cylinder
    c3 = InsideCylinder(center=[1.0, -2.0, 0.5], axis=[1.0, 1.0, 1.0], radius=0.7, length=3.0)
    for x in (SVector(1.2, -1.7, 1.0), SVector(0.0, 0.0, 0.0), SVector(3.0, 0.0, 3.0))
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c3, x), x) ≈ Packmol.constraint_gradient(c3, x)
    end
end
