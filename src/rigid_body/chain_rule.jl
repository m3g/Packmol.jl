#=

This functions apply the chain rule to compute the gradient of the minimum-distance
function relative to the rotations and translation of the rigid-body molecules, from 
the gradient computed from the displacement of the atoms in cartesian coordinates. 

This was implemented originally in the Fortran code of Packmol, by J. M. Martínez.

Extended here for the 2D case.

=#
function gradient_chain_rule!(g::Vector{T}, x, system::PackmolSystem{N}) where {T,N}
    (; ntype, nmols, natoms, comptype, ntotmol) = system
    fill!(zero(T), g)
    k1 = 0
    k2 = ntotmol * N
    icart = 0
    for itype in 1:ntype
        if !comptype(itype)
            icart = icart + nmols(itype) * natoms(itype)
            continue
        end
        for _ in 1:nmols(itype)
            xmol = SVector{N}(x[k1+l] for l in 1:N)
            partial_derivatives!(g, xmol, system)
            k1 += N
            k2 += N
        end
    end
    return g
end

#=

3D case

=#
function partial_derivatives!(g, x::SVector{3}, system::PackmolSystem{3})
    (; gxcar, coor, natoms, idfirst) = system
    beta = x[1]
    gama = x[2]
    teta = x[3]
    sb, cb = sincos(beta)
    sg, cg = sincos(gama)
    st, ct = sincos(teta)
    dv1beta = SVector(-cb * sg * ct - sb * cg, -cb * cg * ct + sb * sg, cb * st)
    dv2beta = SVector(-sb * sg * ct + cb * cg, -sb * cg * ct - cb * sg, sb * st)
    dv3beta = zero(SVector{3,T})
    dv1gama = SVector(-sb * cg * ct - cb * sg, sb * sg * ct - cb * cg, zero(T))
    dv2gama = SVector(cb * cg * ct - sb * sg, -sg * cb * ct - cg * sb, zero(T))
    dv3gama = SVector(cg * st, -sg * st, zero(T))
    dv1teta = SVector(sb * sg * st, sb * cg * st, sb * ct)
    dv2teta = SVector(-cb * sg * st, -cb * cg * st, -cb * ct)
    dv3teta = SVector(sg * ct, cg * ct, -st)
    idatom = idfirst(itype) - 1
    for _ in 1:natoms(itype)
        icart = icart + 1
        idatom = idatom + 1
        for k in 1:3
            g[k1+k] += gxcar[icart][k]
        end
        for k in 1:3
            g[k2+1] += (coor[idatom][1] * dv1beta[k] + coor[idatom][1] * dv2beta[k] + coor[idatom][1] * dv3beta[k]) * gxcar[icart][k]
            g[k2+2] += (coor[idatom][2] * dv1gama[k] + coor[idatom][2] * dv2gama[k] + coor[idatom][2] * dv3gama[k]) * gxcar[icart][k]
            g[k2+3] += (coor[idatom][3] * dv1teta[k] + coor[idatom][3] * dv2teta[k] + coor[idatom][3] * dv3teta[k]) * gxcar[icart][k]
        end
    end
    return g
end

#=

2D case

=#
function partial_derivatives!(g, x::SVector{2}, system::PackmolSystem{2})
    error("not implemented.")
    (; gxcar, coor, natoms, idfirst) = system
    beta = x[1]
    gama = x[2]
    sb, cb = sincos(beta)
    sg, cg = sincos(gama)
    # The following vector definitions are not correct, they are just placeholders
    dv1beta = SVector(-cb * sg * ct - sb * cg, -cb * cg * ct + sb * sg, cb * st)
    dv2beta = SVector(-sb * sg * ct + cb * cg, -sb * cg * ct - cb * sg, sb * st)
    dv1gama = SVector(-sb * cg * ct - cb * sg, sb * sg * ct - cb * cg, zero(T))
    dv2gama = SVector(cb * cg * ct - sb * sg, -sg * cb * ct - cg * sb, zero(T))
    idatom = idfirst(itype) - 1
    for _ in 1:natoms(itype)
        icart = icart + 1
        idatom = idatom + 1
        for k in 1:2
            g[k1+k] += gxcar[icart][k]
        end
        for k in 1:2
            g[k2+1] += (coor[idatom][1] * dv1beta[k] + coor[idatom][1] * dv2beta[k]) * gxcar[icart][k]
            g[k2+2] += (coor[idatom][2] * dv1gama[k] + coor[idatom][2] * dv2gama[k]) * gxcar[icart][k]
        end
    end
    return g
end

@testitem "gradient chain rule 3D" begin

end

@testitem "gradient chain rule 2D" begin

end
