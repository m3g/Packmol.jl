@testitem "run_packmol_legacy" begin
    using Packmol
    using Packmol_jll
    test_dir = Packmol.src_dir*"/../test/input_files"
    run_packmol_legacy("$test_dir/water_box.inp")
    @test isfile("$test_dir/water_box.pdb")
    run_packmol_legacy("$test_dir/ieee_signaling.inp")
    @test isfile("$test_dir/ieee_signaling_box.pdb")
    cd(test_dir)
    rm("water_box.pdb")
    run_packmol_legacy("water_box.inp")
    @test isfile("water_box.pdb")
end

@testitem "run_packmol_legacy without Packmol_jll" begin
    # Run in a fresh process: in the test session Packmol_jll may already be
    # loaded (and the extension active) by another test item.
    project = dirname(Base.active_project())
    code = """
        using Packmol
        try
            run_packmol_legacy("input.inp")
        catch e
            print(e isa ArgumentError && occursin("Packmol_jll", e.msg))
        end
    """
    @test readchomp(`$(Base.julia_cmd()) --project=$project -e $code`) == "true"
end
