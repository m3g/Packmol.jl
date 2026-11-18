# A polymer for the tests: a chain of beads 1.54 Å apart (closer than the tolerance, such
# that consecutive beads of different monomers are excluded pairs), in monomers of 3
# beads, written to a PDB file with one model per conformation: the first all-trans,
# the others with random dihedrals (no beads more than two bonds apart closer than 2.5 Å)
@testmodule PolymerFiles begin
    using LinearAlgebra: norm, cross, normalize
    using Random: Xoshiro
    const dir = mktempdir()
    function place(a, b, c, bond, angle, dihedral)
        bc = normalize(c - b)
        n = normalize(cross(b - a, bc))
        m = cross(n, bc)
        θ, φ = deg2rad(angle), deg2rad(dihedral)
        return c + bond * (-cos(θ) * bc + sin(θ) * cos(φ) * m + sin(θ) * sin(φ) * n)
    end
    function chain(rng, n, trans::Bool)
        x = [[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.11, 1.43, 0.0]]
        while length(x) < n
            y = place(x[end-2], x[end-1], x[end], 1.54, 112.0, trans ? 180.0 : rand(rng, (180.0, 60.0, -60.0)) + 10 * randn(rng))
            if all(norm(y - x[i]) > 2.5 for i in 1:length(x)-2)
                push!(x, y)
            elseif !trans && rand(rng) < 0.1
                x = x[1:max(3, length(x) - 3)]
            end
        end
        return x
    end
    # `nmon` monomers of 3 beads, `nconf` conformations; returns the file
    function polymer_file(name; nmon = 10, nconf = 1, seed = 1)
        rng = Xoshiro(seed)
        file = joinpath(dir, name)
        open(file, "w") do io
            for model in 1:nconf
                nconf > 1 && println(io, "MODEL     ", model)
                x = chain(rng, 3 * nmon, model == 1)
                for (i, xi) in enumerate(x)
                    println(io, rpad("ATOM  " * lpad(i, 5) * "  C   PLY A" * lpad((i - 1) ÷ 3 + 1, 4), 30),
                        join((lpad(string(round(v, digits=3)), 8) for v in xi)), "  1.00  0.00")
                end
                nconf > 1 && println(io, "ENDMDL")
            end
            println(io, "END")
        end
        return file
    end
    # The structure block of `number` polymers of the file
    polymer_block(file; number, extra = "", cyclic = false) = """
        structure $(basename(file))
          number $number
          $extra
          repeat $(count(l -> startswith(l, "ATOM"), eachline(file)) ÷ 3 ÷ max(1, count(l -> startswith(l, "MODEL"), eachline(file))))
            monomer 3
              link 3 1
          end repeat
          $(cyclic ? "cyclic" : "")
        end structure
        """
    # The bonds of each chain of the output (consecutive beads), and the closest pair of
    # beads of different chains, or more than one bead apart in the same chain
    function check_chains(atoms, natoms_chain)
        x = [[a.x, a.y, a.z] for a in atoms]
        bonds = Float64[]
        dmin = Inf
        for i in eachindex(x), j in i+1:length(x)
            ci, cj = (i - 1) ÷ natoms_chain, (j - 1) ÷ natoms_chain
            d = norm(x[i] - x[j])
            if ci == cj && j == i + 1
                push!(bonds, d)
            elseif ci != cj || j > i + 1
                dmin = min(dmin, d)
            end
        end
        return extrema(bonds), dmin
    end
end
