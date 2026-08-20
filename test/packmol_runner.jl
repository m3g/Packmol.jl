@testitem "run_packmol" begin
    using Packmol
    test_dir = Packmol.src_dir*"/../test/input_files"
    run_packmol("$test_dir/water_box.inp")
    @test isfile("$test_dir/water_box.pdb")
    run_packmol("$test_dir/ieee_signaling.inp")
    @test isfile("$test_dir/ieee_signaling_box.pdb")
    cd(test_dir)
    rm("water_box.pdb")
    run_packmol("water_box.inp")
    @test isfile("water_box.pdb")
end
