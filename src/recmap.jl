abstract type RecombinationMap end

haldane(d) = 0.5*(1-exp(-2d))     # distance -> recombination rate
invhaldane(r) = -0.5*log(1 - 2r)  # recombination rate -> distance

struct LinearMap{T} <: RecombinationMap
    r :: T   # per bp recombination rate
end

n_crossovers(rng::AbstractRNG, m::LinearMap, x) = 
    m.r == 0.0 ? 0 : min(length(x), rand(rng, Poisson((length(x)-1) * m.r)))

crossovers(rng::AbstractRNG, _::LinearMap, x, n) = 
    sort(sample(rng, 1:length(x), n, replace=false))

function rand_breakpoints(rng, recmap, tgt)
    n = n_crossovers(rng, recmap, tgt)
    n == 0 && return Int64[]
    breakpoints = crossovers(rng, recmap, tgt, n)
    return breakpoints
end





