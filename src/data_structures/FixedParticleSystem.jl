#
# Cell count for a CellListMap box scales as volume / cutoff^D: for a large,
# sparse box (extent set by where particles *could* be, not by how many there
# are) a naively small cutoff can imply a cell count many orders of magnitude
# larger than the particle count, which is what actually blows up memory —
# the cell array itself has to be allocated at that size before any particles
# are even placed, regardless of parallel/serial construction. This grows
# `min_cutoff` (fewer, bigger cells) just enough to keep the cell count at or
# below `max_particles`, and never shrinks it below `min_cutoff`. Growing the
# cutoff only ever finds more candidate pairs, never fewer, so this is safe
# for correctness as long as callers still filter pairs by the true tolerance
# themselves (as InteratomicDistanceFG and overlaps_fixed do).
function _capped_cutoff(volume::T, min_cutoff::T, max_particles::Real, D::Int) where {T}
    ncells = min(volume / min_cutoff^D, T(max_particles))
    return max(min_cutoff, (volume / ncells)^(one(T) / D))
end

#
# Build a ParticleSystem2 from pre-collected fixed atom positions.
# ypositions = fixed atoms (cell list built once and never rebuilt).
# xpositions = placeholder (first fixed atom); overwritten before each pairwise! call.
# Returns `nothing` if positions are empty.

#
# The output stored in in the ParticleSystem is a compound structure which will 
# be used to compute different properties on each call to pairwise!. The `compute` field of 
# OutputProperty defines if the property will be computed on that call, or not.
#
# Special copy, reset, and reducer functions are defined for each property and invoqued 
# accordingly to the required property computation.
#
mutable struct OutputProperty{T}
    value::T
    compute::Bool
end

# In the case of some special value type, define the appropriate _copy function.
_copy(x) = copy(x)
Base.copy(x::OutputProperty) = OutputProperty(_copy(x.value), x.compute)

@kwdef mutable struct FixedParticleSystemOutput{T}
    overlaps_fixed::OutputProperty{Bool} = OutputProperty(false, false)
    molecule_badness::OutputProperty{Vector{T}} = OutputProperty(T[], false) # output for molecule badness
end

# Set the compute flag to false for all properties
function reset_compute!(x::FixedParticleSystemOutput) 
    for name in propertynames(x)
        setfield!(getfield(x, name), :compute, false)
    end
    return x
end

# Set the compute flag for a property to true
compute_property!(x::FixedParticleSystemOutput, s) = setfield!(getfield(x, s), :compute, true)

CellListMap.copy_output(x::FixedParticleSystemOutput) =
    FixedParticleSystemOutput((_copy(getfield(x, name)) for name in propertynames(x))...)

_reset!(::Val{:overlaps_fixed}, x) = (x.overlaps_fixed.compute) && (x.overlaps_fixed.value = false)
_reset!(::Val{:molecule_badness}, x) = (x.molecule_badness.compute) && (fill!(x.molecule_badness.value, zero(eltype(x.molecule_badness.value))))
function CellListMap.reset_output!(x::FixedParticleSystemOutput)
    for name in propertynames(x)
        _reset!(Val(name), x)
    end
    return x
end 

_reducer!(::Val{:overlaps_fixed}, x, y) = 
    (x.overlaps_fixed.compute) && (x.overlaps_fixed.value = x.overlaps_fixed.value | y.overlaps_fixed.value)
_reducer!(::Val{:molecule_badness}, x, y) = 
    (x.molecule_badness.compute) && (x.molecule_badness.value .+= y.molecule_badness.value)
function CellListMap.reducer!(x::F, y::F) where {F<:FixedParticleSystemOutput}
    for name in propertynames(x)
        _reducer!(Val(name), x, y)
    end
    return x
end

function _build_fixed_particle_system(
    fixed_atom_positions::Vector{SVector{D,T}}, tol::T, nmols::Int;
    nthreads=Threads.nthreads(), parallel::Bool=nthreads > 1,
    unitcell::Union{Nothing,AbstractVecOrMat}=nothing,
) where {D,T}
    isempty(fixed_atom_positions) && return nothing
    output = FixedParticleSystemOutput{T}()
    output.molecule_badness.value = zeros(T, nmols)
    compute_property!(output, :molecule_badness)
    return ParticleSystem(
        xpositions=fixed_atom_positions[1:1],  # placeholder; overwritten before each pairwise! call
        ypositions=fixed_atom_positions,
        cutoff=tol,
        output=output,
        parallel=parallel,
        unitcell=unitcell,
    )
end

#
# Build a ParticleSystem2 from fixed atom positions for efficient overlap checks.
# Returns `nothing` if there are no fixed atoms or avoid_overlap is false.
#
# `nmols` sizes the per-molecule `molecule_badness` output buffer: callers that
# index a whole system's worth of molecules (e.g. molecule_badness.jl, which
# uses a single shared instance) need `packmol_system.nmols`; callers that
# deep-copy one instance per thread/task to check a single trial molecule at a
# time (overlaps_fixed.jl) only ever touch index 1, so pass `nmols=1` there —
# CellListMap replicates the output buffer once per internal batch, and again
# for every deep copy, so an oversized buffer here is multiplied by both
# factors and can blow up memory on large systems.
# `parallel` similarly should be `false` for per-task deep copies (they already
# run inside a caller-managed thread), to avoid CellListMap allocating its own
# internal per-batch output copies on top of that.
#
# `unitcell`, if given, is passed straight through to `ParticleSystem`. Without
# it, CellListMap treats the system as non-periodic and auto-fits a box to
# `limits(xpositions, ypositions)` on every `update!`/`pairwise!` call whose
# trial point falls outside the current box — an O(nmols_fixed) scan, and
# potentially a full cell-list rebuild, repeated on every single trial for
# callers like overlaps_fixed.jl. Passing an explicit box comfortably covering
# every point that will ever be queried makes CellListMap treat it as
# periodic instead, skipping that re-fit entirely.
#
function build_fixed_particle_system(
    packmol_system::PackmolSystem{D,T};
    nmols::Int=packmol_system.nmols, parallel::Bool=Threads.nthreads() > 1,
    unitcell::Union{Nothing,AbstractVecOrMat}=nothing,
    cutoff::T=packmol_system.tolerance,
) where {D,T}
    fixed_atom_positions = _collect_fixed_positions(packmol_system)
    return _build_fixed_particle_system(fixed_atom_positions, cutoff, nmols; parallel, unitcell)
end

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
