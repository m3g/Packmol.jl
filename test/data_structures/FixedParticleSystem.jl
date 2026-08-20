@testitem "FixedParticleSystem" begin
    using Packmol
    using Packmol: FixedParticleSystemOutput, 
                   OutputProperty,
                   compute_property!,
                   reset_compute!
    using CellListMap: copy_output,
                       reset_output!,
                       reducer!

    T = Float64

    #
    # Test setup of particle system
    #
    o = FixedParticleSystemOutput{T}()
    @test o.overlaps_fixed.compute == false
    @test o.molecule_badness.compute == false
    @test o.molecule_badness.value == T[]

    compute_property!(o, :molecule_badness)
    @test o.overlaps_fixed.compute == false
    @test o.molecule_badness.compute == true
    @test o.molecule_badness.value == T[]

    reset_compute!(o)
    @test o.overlaps_fixed.compute == false
    @test o.molecule_badness.compute == false

    #
    # Test parallel interface functions
    #
    o1 = FixedParticleSystemOutput{T}()
    o2 = FixedParticleSystemOutput{T}()
    compute_property!(o1, :overlaps_fixed)
    compute_property!(o1, :molecule_badness)

    o1.overlaps_fixed.value = false
    o2.overlaps_fixed.value = true
    o1.molecule_badness.value = [1.0, 2.0]
    o2.molecule_badness.value = [1.0, 2.0]

    reducer!(o1, o2)
    @test o1.overlaps_fixed.value == true
    @test o1.molecule_badness.value == [2.0, 4.0]
    o3 = copy_output(o1)
    @test o3.overlaps_fixed.value == true
    @test o3.molecule_badness.value == [2.0, 4.0]
    reset_output!(o3)
    @test o3.overlaps_fixed.value == false
    @test o3.molecule_badness.value == [0.0, 0.0]

    #
    # Test updating only of the properties
    #
    reset_compute!(o1)
    compute_property!(o1, :molecule_badness)
    o1.overlaps_fixed.value = false
    o2.overlaps_fixed.value = true
    o1.molecule_badness.value = [1.0, 2.0]
    o2.molecule_badness.value = [1.0, 2.0]
    reducer!(o1, o2)
    @test o1.overlaps_fixed.value == false
    @test o1.molecule_badness.value == [2.0, 4.0]
    o3 = copy_output(o1)
    @test o3.overlaps_fixed.value == false
    @test o3.molecule_badness.value == [2.0, 4.0]
    o3.overlaps_fixed.value = true
    reset_output!(o3)
    @test o3.overlaps_fixed.value == true
    @test o3.molecule_badness.value == [0.0, 0.0]

end
