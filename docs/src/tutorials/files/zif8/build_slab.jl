#
# A slab of ZIF-8 (zif8.pdb): 3 x 3 unit cells in x and y (periodic), 2 unit cells thick
# in z, with whole linkers at its faces, from the unit cell of COD entry 7111973
# (Shekhah et al., Chem. Commun. 50, 2089 (2014)) filled by Open Babel:
#
#   obabel 7111973.cif -O zif8_uc.pdb --fillUC strict
#   julia build_slab.jl zif8_uc.pdb zif8.pdb
#
using PDBTools, LinearAlgebra

# Minimum-image distance, with the periods L (Inf: not periodic)
wrap(d, L) = isinf(L) ? d : d - L * round(d / L)
dist(p, q, L) = norm([wrap(p[i] - q[i], L[i]) for i in 1:3])

function zif8_slab(ucfile, outfile; a = 17.033, n = 3)
    uc = read_pdb(ucfile)
    x = [collect(Float64, position(at)) for at in uc]
    el = [element(at) for at in uc]

    # The hydrogens of the methyl groups are disordered in the crystal: two sets of 3,
    # half occupied. One set per carbon is kept: three hydrogens 1.57 Å apart (C–H bonds
    # of 0.96 Å in the CIF); those of the other set are 0.9 and 1.8 Å from them.
    L = (a, a, a)
    remove = Int[]
    for i in eachindex(x)
        el[i] == "C" || continue
        hs = [j for j in eachindex(x) if el[j] == "H" && dist(x[i], x[j], L) < 1.2]
        length(hs) == 6 || continue
        keep = [hs[1]; [h for h in hs[2:end] if 1.4 < dist(x[h], x[hs[1]], L) < 1.7]]
        length(keep) == 3 || error("methyl $i: no set of 3 hydrogens")
        append!(remove, setdiff(hs, keep))
    end
    keepat = setdiff(eachindex(x), remove)
    x = x[keepat]
    el = el[keepat]
    println("unit cell: ", length(x), " atoms") # Zn12 (C4H5N2)24: 276 atoms

    # n x n cells in x and y (periodic), 3 cells stacked in z (not periodic)
    X = Vector{Float64}[]
    E = String[]
    for i in 0:n-1, j in 0:n-1, k in 0:2, (p, e) in zip(x, el)
        push!(X, p + [i * a, j * a, k * a])
        push!(E, e)
    end
    Lxy = (n * a, n * a, Inf)

    # The linkers: the groups of bonded atoms (bonds shorter than 1.75 Å) other than zinc
    nat = length(X)
    linker = zeros(Int, nat)
    nlinkers = 0
    for i in 1:nat
        (E[i] == "Zn" || linker[i] != 0) && continue
        nlinkers += 1
        linker[i] = nlinkers
        stack = [i]
        while !isempty(stack)
            k = pop!(stack)
            for j in 1:nat
                if E[j] != "Zn" && linker[j] == 0 && dist(X[k], X[j], Lxy) < 1.75
                    linker[j] = nlinkers
                    push!(stack, j)
                end
            end
        end
    end

    # The zinc atoms of the two cells in the middle of the stack, and the linkers bound to
    # them (those at the faces are bound to a single zinc)
    zn = [i for i in 1:nat if E[i] == "Zn" && a <= X[i][3] < 3a]
    bound = Set(linker[j] for i in zn for j in 1:nat if E[j] == "N" && dist(X[i], X[j], Lxy) < 2.3)
    slab = sort([zn; [j for j in 1:nat if linker[j] in bound]])

    # The slab from z = 0.5, with x and y wrapped into the periodic cell
    z0 = minimum(X[i][3] for i in slab)
    atoms = Atom{Nothing}[]
    for i in slab
        p = X[i] - [0.0, 0.0, z0 - 0.5]
        push!(atoms, Atom(; name = E[i], resname = "ZIF", resnum = 1, chain = "A",
            x = mod(p[1], n * a), y = mod(p[2], n * a), z = p[3], pdb_element = E[i]))
    end
    println("slab: ", length(atoms), " atoms, ", length(zn), " Zn, ", length(bound),
        " linkers, z from 0.5 to ", round(maximum(at.z for at in atoms); digits = 2))
    write_pdb(outfile, atoms)
end

zif8_slab(ARGS[1], ARGS[2])
