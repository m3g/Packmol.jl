module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import Packmol: test_file, tutorial_file

test_file(name::AbstractString) = joinpath(artifact"packmol_test_files", name)
tutorial_file(name::AbstractString) = joinpath(artifact"packmol_tutorial_files", name)

end
