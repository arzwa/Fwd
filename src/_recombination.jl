abstract type RecombinationMap end

# Haldane's mapping function
# distance -> recombination rate
recrate(d) = 0.5*(1-exp(-2d))     

# recombination rate -> distance
distance(r) = -0.5*log(1 - 2r)  

# recombination rate matrix
rec_matrix(x) = [recrate(abs(x[i] - x[j])) for i=1:length(x), j=1:length(x)]

struct LinearMap{T} <: RecombinationMap
    r :: T   # per bp recombination rate
end

function rand_breakpoints(rng, recmap, tgt)
    n = n_crossovers(rng, recmap, tgt)
    n == 0 && return Int64[]
    breakpoints = crossovers(rng, recmap, tgt, n)
    return breakpoints
end

function n_crossovers(rng::AbstractRNG, m::LinearMap, x) 
    m.r == 0.0 && return 0
    return min(length(x), rand(rng, Poisson((length(x)-1) * m.r)))
end

function crossovers(rng::AbstractRNG, _::LinearMap, x, n)
    sort(sample(rng, 1:length(x), n, replace=false))
end

