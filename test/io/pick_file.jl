@testitem "pick_file" begin
    using Packmol: pick_file, _file_dialog_cmd, _popen_read, _powershell_encode
    withenv("PACKMOL_GUI" => "false") do
        @test_throws ArgumentError pick_file()
        @test_throws ArgumentError packmol()
    end
    cmd = _file_dialog_cmd()
    if Sys.isapple() || Sys.iswindows() || !isnothing(Sys.which("zenity")) || !isnothing(Sys.which("kdialog"))
        @test cmd isa String
    else
        @test isnothing(cmd)
    end
    @test strip(_popen_read("echo packmol")) == "packmol"
    @test _popen_read(Sys.iswindows() ? "exit 1" : "false") == ""
    # Base64 of the UTF-16LE bytes, as PowerShell's -EncodedCommand expects
    @test _powershell_encode("ab") == "YQBiAA=="
    @test _powershell_encode("a") == "YQA="
    @test _powershell_encode("abc") == "YQBiAGMA"
    @test _powershell_encode("\$dialog = é") == "JABkAGkAYQBsAG8AZwAgAD0AIADpAA=="
end

@testitem "terminal output" begin
    using Packmol: _print, _println, @_printf, _use_core_streams
    file = tempname()
    open(file, "w") do io
        redirect_stdout(io) do
            _print("a", 1, ' ')
            _println(2.5, :b)
            @_printf("%5.2f|%s\n", 1.0, "x")
        end
    end
    @test read(file, String) == "a1 2.5b\n 1.00|x\n"
    @test _use_core_streams[] == false
end
