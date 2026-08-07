#
# Lazy, locally-compiled build of the GENCAN FFI shared library.
#
# This is an experimental, local-build-only integration (see src/gencan/
# for the design rationale): the vendored Fortran sources are compiled with
# gfortran on first use and cached in a Scratch.jl directory. There is no
# BinaryBuilder/JLL artifact yet — using `optimizer=:gencan` requires a
# working `gfortran` on the machine running Packmol.jl.
#

const GENCAN_FORTRAN_DIR = joinpath(@__DIR__, "fortran")
const GENCAN_SOURCES = (
    joinpath(GENCAN_FORTRAN_DIR, "callback_bridge.f90"),
    joinpath(GENCAN_FORTRAN_DIR, "gencan_ieee_signal_routines.f90"),
    joinpath(GENCAN_FORTRAN_DIR, "gencan_core.f"),
    joinpath(GENCAN_FORTRAN_DIR, "ahestetic.f90"),
)

const _gencan_build_lock = ReentrantLock()
_gencan_libhandle = Ref{Union{Nothing,Ptr{Cvoid}}}(nothing)

gencan_gfortran_available() = Sys.which("gfortran") !== nothing

#
# Ensure the GENCAN shared library is built and loaded, returning a dlopen
# handle. Rebuilds if missing or if any source file is newer than the
# cached binary. `ccall`'s library reference must be a compile-time
# constant, which a lazily-built path is not — so we dlopen once and
# resolve symbols via Libdl.dlsym (see ffi.jl), the standard pattern for
# dynamically-located libraries.
#
function ensure_libgencan_built()
    lock(_gencan_build_lock) do
        if !gencan_gfortran_available()
            error(
                "optimizer=:gencan requires a `gfortran` compiler on PATH " *
                "(local-build-only integration, no precompiled binary yet). " *
                "Install gfortran or use the default optimizer=:spgbox.",
            )
        end

        scratch_dir = Scratch.@get_scratch!("gencan_build")
        libpath = joinpath(scratch_dir, "libgencan." * Libdl.dlext)

        needs_build = !isfile(libpath) || any(
            src -> mtime(src) > mtime(libpath), GENCAN_SOURCES
        )

        if needs_build
            _gencan_libhandle[] !== nothing && Libdl.dlclose(_gencan_libhandle[])
            _gencan_libhandle[] = nothing
            modcache = mktempdir()
            cmd = Cmd([
                "gfortran", "-shared", "-fPIC", "-O2", "-o", libpath,
                GENCAN_SOURCES..., "-J", modcache,
            ])
            run(cmd)
        end

        if _gencan_libhandle[] === nothing
            _gencan_libhandle[] = Libdl.dlopen(libpath)
        end

        return _gencan_libhandle[]
    end
end

function gencan_libhandle()
    handle = _gencan_libhandle[]
    return handle === nothing ? ensure_libgencan_built() : handle
end

function gencan_symbol(name::Symbol)
    return Libdl.dlsym(gencan_libhandle(), name)
end
