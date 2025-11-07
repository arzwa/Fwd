# NOTE: would it be nicer to work with pairwise recombination rates, as
# that accounts for unlinked loci automatically? But not sure how to sample
# breakpoints efficiently in general.

abstract type RecombinationMap end

struct LinearMap{T} <: RecombinationMap
    maplength :: T  # maplength in Morgans, i.e. expected # of crossovers
end
Base.length(m::LinearMap) = m.maplength

@with_kw struct LinearPhysMap{T} <: RecombinationMap
    physlength :: Int
    maplength  :: T
    rbp        :: T = maplength/physlength   # M/bp
end
Base.length(m::LinearPhysMap) = m.physlength
recrate(m::LinearPhysMap, x, y) = recrate(m.rbp * abs(x - y))

struct Unlinked <: RecombinationMap end
Base.length(m::Unlinked) = 1

"""
    Chromosomes

A genetic map consisting of multiple chromosomes. We assume that loci are
given map positions along the genome as if chromosomes were lined up one
after another. For instance, if we have a genetic map defined by
```
M = Chromosomes([LinearMap(0.7), LinearMap(0.3)])
```
A locus at map position 0.1 on the second chromosome should have coordinate 
in an `Architecture` object `x=0.8`.

To model a set of `L` unlinked loci, use
```
L = 10
M = Chromosomes([Unlinked() for _=1:L])
```
This is equivalent to 
```
M = Chromosomes([LinearPhysMap(physlength=1, maplength=0.0) for _=1:L])
```
The latter however admits mixing linked with unlinked stuff.
"""
struct Chromosomes{M<:RecombinationMap} <: RecombinationMap 
    maps :: Vector{M}
end    

maplength(m) = m.maplength

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

function rand_breakpoints(rng, m::LinearPhysMap)
    n = rand(rng, Poisson(maplength(m)))
    bps = sample(rng, 1:m.physlength-1, n, replace=false)
    [sort!(bps); m.physlength]
end

rand_breakpoints(_, m::Unlinked) = [1]

# XXX could have a specialized implementation for unlinked architectures
function rand_breakpoints(rng, m::Chromosomes)
    C = 0.0
    bps = map(m.maps) do recmap
        bps = rand_breakpoints(rng, recmap)
        bps .+= C
        C += length(recmap)
        bps
    end
    for chrom in bps[1:end-1]
        rand(rng) < 0.5 && pop!(chrom)
        # the last entry for each chromosome is the chromosome endpoint,
        # if we keep it among breakpoints, there's a recombination between
        # unlinked chromosomes, if we remove it, there is no recombination.
    end
    vcat(bps...)
end

# `recombine!` is a general function, different sorts of genetic map should
# implement their specific `rand_breakpoints` function. 
"""
    recombine!(z, breakpoints, x, y, xs)

Recombine `x` and `y` assuming crossover recombination at `breakpoints`,
assuming the entries of `x` and `y` are at map positions `xs` (should be sorted),
write to `z`.

!!! note: This function is deterministic, for a given set of breakpoints and
`x` and `y` haplotypes, it will always return the same recombinant haplotype.
To obtain a random recombinant haplotype for a given set of brekapoints and
haplotype (i.e. a random pick of the two recombinant haplotypes), one should
randomize the order of the `x` and `y` arguments.
"""
function recombine!(z, breakpoints, x, y, xs, onx=true) 
    length(z) == 0 && return
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

