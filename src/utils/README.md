# utils

Basic functions used by the rest of the code. All are compatible with `juliac --trim`
(see the notes in `CLAUDE.md`): in a trimmed executable `stdout` and `stderr` are not
usable through the usual functions, and `@sync` is not compiled.

| File | Contents |
|---|---|
| `terminal.jl` | Printing to the terminal: `_print`, `_println`, `@_printf`, `_with_stdout`/`_with_stderr`, `_readline`. In the executable they write to the `Core` streams. |
| `progress_bar.jl` | A minimal progress bar (`ProgressBar`, `next!`, `finish!`), replacing ProgressMeter.jl. |
| `threads.jl` | `@_sync`, a trim-safe `@sync`, and `_ntasks`, the number of tasks of the parallel loops. |
| `pick_file.jl` | `pick_file`: a file dialog with the tools of the operating system (osascript, PowerShell, zenity, kdialog), run with libc `popen`, so that no GUI toolkit is bundled. |
