"""
    pick_file()

Opens a file dialog to choose the Packmol input file, and returns its path.

The dialog is provided by tools that come with the operating system — `osascript`
on macOS, PowerShell on Windows, and `zenity` or `kdialog` on Linux — so no GUI
toolkit is bundled with Packmol. If none of them is available, the path is read
from the terminal instead.

To disable the file dialog, set the environment variable `PACKMOL_GUI="false"`.

"""
function pick_file()
    if get(ENV, "PACKMOL_GUI", "true") == "false"
        throw(ArgumentError("""\n
            Environment variable PACKMOL_GUI is set to false. Provide the input file explicitly,
            or set PACKMOL_GUI to true to use the file dialog.

        """))
    end
    cmd = _file_dialog_cmd()
    if isnothing(cmd)
        _print("Input file: ")
        input_file = strip(_readline())
    else
        # Cancelling the dialog makes the tool exit with a nonzero status.
        input_file = strip(_popen_read(cmd))
    end
    isempty(input_file) && throw(ArgumentError("No input file selected."))
    return String(input_file)
end

const _dialog_title = "Select the Packmol input file"

# The shell command that opens the file dialog, or `nothing` if there is no tool
# for it. A command line string, rather than a `Cmd`, because it is run with
# `_popen_read` (see there).
function _file_dialog_cmd()
    if Sys.isapple()
        return "osascript -e 'POSIX path of (choose file with prompt \"$_dialog_title\")'"
    elseif Sys.iswindows()
        # The hidden TopMost form keeps the dialog from opening behind the terminal.
        script = join((
            "[Console]::OutputEncoding = [System.Text.Encoding]::UTF8",
            "Add-Type -AssemblyName System.Windows.Forms",
            "\$owner = New-Object System.Windows.Forms.Form -Property @{TopMost = \$true}",
            "\$dialog = New-Object System.Windows.Forms.OpenFileDialog",
            "\$dialog.Title = '$_dialog_title'",
            "\$dialog.Filter = 'Packmol input (*.inp)|*.inp|All files (*.*)|*.*'",
            "if (\$dialog.ShowDialog(\$owner) -eq 'OK') { \$dialog.FileName }",
        ), "; ")
        # Passed encoded (base64 of UTF-16LE), so that it is not parsed by the shell.
        return "powershell -NoProfile -STA -EncodedCommand $(_powershell_encode(script))"
    elseif !isnothing(Sys.which("zenity"))
        return "zenity --file-selection --title='$_dialog_title' --file-filter='*.inp' --file-filter='*'"
    elseif !isnothing(Sys.which("kdialog"))
        return "kdialog --title '$_dialog_title' --getopenfilename . '*.inp'"
    end
    return nothing
end

#
# Run the shell command `command`, and return what it writes to stdout, or an empty
# string if it fails. Uses the C library's `popen`, because Julia's own process
# spawning (`run`, `read(::Cmd)`) does not compile with `juliac --trim`.
#
function _popen_read(command::String)
    fp = @static if Sys.iswindows()
        ccall(:_popen, Ptr{Cvoid}, (Cstring, Cstring), command, "r")
    else
        ccall(:popen, Ptr{Cvoid}, (Cstring, Cstring), command, "r")
    end
    fp == C_NULL && return ""
    output = IOBuffer()
    chunk = Vector{UInt8}(undef, 4096)
    while true
        n = ccall(:fread, Csize_t, (Ptr{UInt8}, Csize_t, Csize_t, Ptr{Cvoid}), chunk, 1, length(chunk), fp)
        n == 0 && break
        write(output, view(chunk, 1:Int(n)))
    end
    status = @static if Sys.iswindows()
        ccall(:_pclose, Cint, (Ptr{Cvoid},), fp)
    else
        ccall(:pclose, Cint, (Ptr{Cvoid},), fp)
    end
    return status == 0 ? String(take!(output)) : ""
end

# The base64 encoding of the UTF-16LE bytes of `script`, as PowerShell's
# `-EncodedCommand` expects.
function _powershell_encode(script::String)
    bytes = UInt8[]
    for u in transcode(UInt16, script)
        push!(bytes, UInt8(u & 0xff), UInt8(u >> 8))
    end
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
    encoded = IOBuffer()
    for i in 1:3:length(bytes)
        nb = min(3, length(bytes) - i + 1)
        b = UInt32(bytes[i]) << 16
        nb > 1 && (b |= UInt32(bytes[i+1]) << 8)
        nb > 2 && (b |= UInt32(bytes[i+2]))
        for j in 1:4
            if j <= nb + 1
                write(encoded, alphabet[((b >> (6 * (4 - j))) & 0x3f) + 1])
            else
                write(encoded, '=')
            end
        end
    end
    return String(take!(encoded))
end
