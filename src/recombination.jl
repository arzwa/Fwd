abstract type RecombinationMap end

# would be nicer to work with pairwise recombination rates, as that accounts
# for unlinked loci automatically, but not sure how to sample breakpoints
# efficiently in general.

struct Unlinked <: RecombinationMap end
maplength(m::Unlinked) = Inf

struct LinearMap{T} <: RecombinationMap
    maplength :: T  # maplength in Morgans, i.e. expected # of crossovers
end

maplength(m::LinearMap) = m.maplength

# Haldane's mapping function
# distance -> recombination rate
recrate(d) = 0.5*(1-exp(-2d))     

# recombination rate -> distance
distance(r) = -0.5*log(1 - 2r)  

# recombination rate matrix
rec_matrix(x) = [recrate(abs(x[i] - x[j])) for i=1:length(x), j=1:length(x)]

function rand_breakpoints(rng, m::LinearMap)
    L = maplength(m)
    n = rand(rng, Poisson(L))
    bps = rand(rng, n) .* L
    [sort!(bps) ; L]
end

# `recombine!` is a general function, different sorts of genetic map should
# implement their specific `rand_breakpoints` function. 
"""
    recombine!(z, breakpoints, x, y, xs)

Recombine `x` and `y` assuming crossover recombination at `breakpoints`,
assuming the entries of `x` and `y` are at map positions `xs`, write to `z`.

!!! note: This function is deterministic, for a given set of breakpoints and
`x` and `y` haplotypes, it will always return the same recombinant haplotype.
To obtain a random recombinant haplotype for a given set of brekapoints and
haplotype (i.e. a random pick of the two recombinant haplotypes), one should
randomize the order of the `x` and `y` arguments.
"""
function recombine!(z, breakpoints, x, y, xs)
    length(z) == 0 && return
    onx = true
    i   = 1
    for bp in breakpoints
        while i <= length(xs) && xs[i] <= bp
            z[i] = onx ? x[i] : y[i]
            i += 1
        end
        i > length(xs) && break
        xs[i] > bp && (onx = !onx)
    end
    i > length(z) && return z
    z[i:end] .= onx ? x[i:end] : y[i:end]
    return z
end

function recombination!(target, rng, recmap, x1, x2, arch)
    breakpoints = rand_breakpoints(rng, recmap)
    recombine!(target, breakpoints, x1, x2, arch.xs)
end

function recombination!(target, rng, recmap::Unlinked, x1, x2, _)
    for i in 1:length(target)
        target[i] = rand(rng) < 0.5 ? x1[i] : x2[i] 
    end
end


