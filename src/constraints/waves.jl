#
# Wave constraints: above/below a sinusoidal surface. `cos` is implemented as
# `sin` with the phase shifted by pi/2 (`cos(t) == sin(t + pi/2)`), so there
# is a single underlying struct for both.
#
# Two flavors, both generalizing `Plane` (they degenerate exactly to it when
# `amplitude == 0`), mirroring `Gaussian`/`RadialGaussian` in gaussians.jl:
#
#   - `Wave` (ridge): the surface is `up·x = h(along·x)`, a sine wave
#     extruded along the direction perpendicular to `along` (within the
#     plane orthogonal to `up`). `up` is used un-normalized, as in `Plane`;
#     `along` is normalized internally.
#   - `RadialWave` (concentric ripples): the surface is `up·x = h(r)`, where
#     `r` is the radial distance from the axis line through `center`
#     parallel to `up`. `up` is normalized internally, since it doubles as
#     the axis for the radial decomposition (as in `Cylinder`). Unlike
#     `RadialGaussian`, `h` is evaluated on `r` itself (not `r^2`), so the
#     gradient has a removable singularity at `r == 0` (the ring's center),
#     handled by treating the radial direction there as zero.
#
export Wave, AboveSin, BelowSin, AboveCos, BelowCos
export RadialWave, AboveRadialSin, BelowRadialSin, AboveRadialCos, BelowRadialCos

# Default weights
weight_default[:wave] = 5.0
weight_default[:radial_wave] = 5.0

#
# Ridge wave
#
@kwdef struct Wave{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    along::SVector{3,T}
    d0::T
    amplitude::T
    wavelength::T
    phase::T
    weight::T = weight_default[:wave]
end

function Wave{Placement}(up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) where {Placement}
    T = promote_type(eltype(up), eltype(along), typeof(d0), typeof(amplitude), typeof(wavelength), typeof(phase), typeof(weight))
    return Wave{Placement,T}(SVector{3,T}(up), SVector{3,T}(along), T(d0), T(amplitude), T(wavelength), T(phase), T(weight))
end
Wave{Placement}(; up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) where {Placement} =
    Wave{Placement}(up, along, d0, amplitude, wavelength, phase, weight)

AboveSin(args...; kargs...) = Wave{Over}(args...; kargs...)
BelowSin(args...; kargs...) = Wave{Below}(args...; kargs...)
AboveCos(up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) =
    Wave{Over}(up, along, d0, amplitude, wavelength, phase + oftype(phase, pi) / 2, weight)
AboveCos(; up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) =
    AboveCos(up, along, d0, amplitude, wavelength, phase, weight)
BelowCos(up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) =
    Wave{Below}(up, along, d0, amplitude, wavelength, phase + oftype(phase, pi) / 2, weight)
BelowCos(; up, along, d0, amplitude, wavelength, phase, weight=weight_default[:wave]) =
    BelowCos(up, along, d0, amplitude, wavelength, phase, weight)

# Profile value and derivative at s = along_unit·x
function _wave_h(c::Wave, s)
    (; d0, amplitude, wavelength, phase) = c
    k = 2 * oftype(s, pi) / wavelength
    return d0 + amplitude * sin(k * s + phase)
end
function _wave_dh(c::Wave, s)
    k = 2 * oftype(s, pi) / c.wavelength
    return c.amplitude * k * cos(k * s + c.phase)
end

function _wave_v_gradv(c::Wave, x)
    a = c.along / norm(c.along)
    s = dot(a, x)
    h = _wave_h(c, s)
    v = dot(c.up, x) - h
    dh = _wave_dh(c, s)
    gradv = c.up - dh * a
    return v, gradv
end

_op_wave(::Wave{Over}, v) = v < zero(v)
_op_wave(::Wave{Below}, v) = v > zero(v)

function constraint_penalty(c::Wave, x)
    v, _ = _wave_v_gradv(c, x)
    if _op_wave(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Wave, x)
    v, gradv = _wave_v_gradv(c, x)
    if _op_wave(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Radial wave (concentric ripples)
#
@kwdef struct RadialWave{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    center::SVector{3,T}
    d0::T
    amplitude::T
    wavelength::T
    phase::T
    weight::T = weight_default[:radial_wave]
end

function RadialWave{Placement}(up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) where {Placement}
    T = promote_type(eltype(up), eltype(center), typeof(d0), typeof(amplitude), typeof(wavelength), typeof(phase), typeof(weight))
    return RadialWave{Placement,T}(SVector{3,T}(up), SVector{3,T}(center), T(d0), T(amplitude), T(wavelength), T(phase), T(weight))
end
RadialWave{Placement}(; up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) where {Placement} =
    RadialWave{Placement}(up, center, d0, amplitude, wavelength, phase, weight)

AboveRadialSin(args...; kargs...) = RadialWave{Over}(args...; kargs...)
BelowRadialSin(args...; kargs...) = RadialWave{Below}(args...; kargs...)
AboveRadialCos(up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) =
    RadialWave{Over}(up, center, d0, amplitude, wavelength, phase + oftype(phase, pi) / 2, weight)
AboveRadialCos(; up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) =
    AboveRadialCos(up, center, d0, amplitude, wavelength, phase, weight)
BelowRadialCos(up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) =
    RadialWave{Below}(up, center, d0, amplitude, wavelength, phase + oftype(phase, pi) / 2, weight)
BelowRadialCos(; up, center, d0, amplitude, wavelength, phase, weight=weight_default[:radial_wave]) =
    BelowRadialCos(up, center, d0, amplitude, wavelength, phase, weight)

# Decompose x - center into the axial coordinate `w` (along the normalized
# `up`) and the perpendicular vector `perp`, whose norm `r` is the radial
# distance to the up-axis line through `center`.
function _radial_wave_axial_perp(c::RadialWave, x)
    u = c.up / norm(c.up)
    a = x - c.center
    w = dot(a, u)
    perp = a - w * u
    return u, w, perp
end

function _radial_wave_v_gradv(c::RadialWave, x)
    u, w, perp = _radial_wave_axial_perp(c, x)
    r = norm(perp)
    k = 2 * oftype(r, pi) / c.wavelength
    h = c.d0 + c.amplitude * sin(k * r + c.phase)
    v = w - h
    dh_dr = c.amplitude * k * cos(k * r + c.phase)
    # Radial direction is undefined exactly at r == 0; treat its
    # contribution to the gradient as zero there (removable singularity).
    rdir = iszero(r) ? zero(perp) : perp / r
    gradv = u - dh_dr * rdir
    return v, gradv
end

_op_wave(::RadialWave{Over}, v) = v < zero(v)
_op_wave(::RadialWave{Below}, v) = v > zero(v)

function constraint_penalty(c::RadialWave, x)
    v, _ = _radial_wave_v_gradv(c, x)
    if _op_wave(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::RadialWave, x)
    v, gradv = _radial_wave_v_gradv(c, x)
    if _op_wave(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["above sin"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above sin' constraint data for $(structure_data[:filename]).")
    end
    return Wave{Over,T}(; up, along, d0, amplitude, wavelength, phase)
end
parse_constraint["below sin"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below sin' constraint data for $(structure_data[:filename]).")
    end
    return Wave{Below,T}(; up, along, d0, amplitude, wavelength, phase)
end

parse_constraint["above cos"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above cos' constraint data for $(structure_data[:filename]).")
    end
    return Wave{Over}(up, along, d0, amplitude, wavelength, phase + T(pi) / 2)
end
parse_constraint["below cos"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below cos' constraint data for $(structure_data[:filename]).")
    end
    return Wave{Below}(up, along, d0, amplitude, wavelength, phase + T(pi) / 2)
end

parse_constraint["above radial_sin"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above radial_sin' constraint data for $(structure_data[:filename]).")
    end
    return RadialWave{Over,T}(; up, center, d0, amplitude, wavelength, phase)
end
parse_constraint["below radial_sin"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below radial_sin' constraint data for $(structure_data[:filename]).")
    end
    return RadialWave{Below,T}(; up, center, d0, amplitude, wavelength, phase)
end

parse_constraint["above radial_cos"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above radial_cos' constraint data for $(structure_data[:filename]).")
    end
    return RadialWave{Over}(up, center, d0, amplitude, wavelength, phase + T(pi) / 2)
end
parse_constraint["below radial_cos"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, wavelength, phase = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below radial_cos' constraint data for $(structure_data[:filename]).")
    end
    return RadialWave{Below}(up, center, d0, amplitude, wavelength, phase + T(pi) / 2)
end
