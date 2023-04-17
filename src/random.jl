import Random
using StableRNGs
function init_random(options)
    if options.seed > 0
        if options.StableRNG == true
            RNG = StableRNGs.StableRNG(options.seed)
        else
            RNG = Random.Xoshiro(options.seed)
        end
    else
        if options.StableRNG == true
            seed = abs(rand(Int))
            RNG = StableRNGs.StableRNG(seed)
        else
            RNG = Random.Xoshiro()
        end
    end
    return RNG
end
random(RNG) = rand(RNG)
random(RNG, arg) = rand(RNG, arg)
