# Test files: the `packmol` set of https://github.com/m3g/MolSimTestFiles (a lazy artifact,
# see Artifacts.toml). The tests use a writable copy, since they write output files next to
# the input files, and the artifact is read-only.
@testmodule TestFiles begin
    using LazyArtifacts: LazyArtifacts
    import Packmol
    const dir = mktempdir()
    cp(Packmol.test_file(""), dir; force=true)
    chmod(dir, 0o755; recursive=true)
    test_file(name::AbstractString) = joinpath(dir, name)
end
