# Architecture is a struct which should define functions giving relevant
# parameters at a site
abstract type Locus end

struct HaploidLocus{T} <: Locus
    s :: T  # selection coefficient
    u :: T  # mutation rate
end

isneutral(l::HaploidLocus) = l.s == 0.0

struct Region{V<:AbstractVector}
    xs :: V              
    L  :: Int64
    ls :: Vector{Int64}
end

Region(x::AbstractVector{Int}) = Region(x, length(x), Int64[])
Region(x::AbstractVector{<:AbstractVector}) = Region(x, sum(length.(x)), length.(x))

## General architecture: Have a set of 'kinds of loci' and regions (which may be
## disjoint) with loci of these types (somewhat like SLiM).
struct Architecture{L<:Locus,R<:Region}
    loci :: Vector{L}   # different types of loci
    xs   :: Vector{R}   # locations/regions of loci
end

# XXX This is not sufficiently sparse!
# we should iterate over the nonzero indices
function fitness(A::Architecture, genome)
    @unpack loci, xs = A
    logw = 0.0 
    for (x, l) in zip(xs, loci)
        logw += logfitness(l, x, genome)
    end
    exp(logw)
end

function logfitness(l::HaploidLocus, x::Region, genome) 
    l.s == 0.0 && return 0.0
    logw = 0.0
    for subregion in x.xs
        #logw += sum(genome[subregion]) * l.s
        logw += sum([genome[r] for r in subregion]) * l.s  # faster
    end
    return logw
end

function mutation!(rng::AbstractRNG, tgt, arch::Architecture)
    @unpack xs, loci = arch
    for (x, l) in zip(xs, loci)
        U  = l.u * x.L  # l is a locus, x is a region
        nm = rand(rng, Poisson(U))
        for _=1:nm
            randmut!(rng, x, tgt)
        end
    end
    dropzeros!(tgt)  # important!
    return tgt
end

# Region is a list of positions or unitrange or similar
function randmut!(rng::AbstractRNG, region::Region{<:AbstractVector{Int}}, tgt)
    i = rand(rng, region.xs)
    tgt[i] = !tgt[i]
end

# Region is a list of unitranges or similar
function randmut!(rng::AbstractRNG, region::Region, tgt)
    i = sample(rng, 1:length(region.xs), Weights(region.ls))
    j = rand(rng, region.xs[i])
    tgt[j] = !tgt[j]
end

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

"""
    recombination(tgt, src1, src2, breakpoints)

Recombine (crossover recombination) haplotypes `src[1]` and `src[2]` and write
to `tgt`, for a given set of breakpoints.

Note that we always start with `src1`, so in applications it is assumed that
`src1` and `src2` are randomized already.
"""
function recombination(src1, src2, breakpoints)
    length(breakpoints) == 0 && return copy(src1) 
    # we walk along the breakpoints and nonzero indices.
    L   = length(src1)
    xs  = (src1.nzind, src2.nzind)
    x3  = empty(xs[1])
    bps = [breakpoints ; L]
    ks  = [1, 1]
    for (i,bp) in enumerate(bps)
        f, b = isodd(i) ? (1, 2) : (2, 1)
        while ks[f] <= length(xs[f]) && xs[f][ks[f]] <= bp
            push!(x3, xs[f][ks[f]])
            ks[f] += 1
        end
        while ks[b] <= length(xs[b]) && xs[b][ks[b]] <= bp
            ks[b] += 1
        end
    end
    n = length(x3)
    SparseVector(L, x3, ones(Bool, n))
end

# ----------------------------------------------------------------

#struct NeutralArchitecture{T} <: Architecture 
#    u :: T   # mutation rate 
#end
#
## This approach allows to model arbitrary mutation rate variation
## (assume everything that is not in `loci` is neutral with neutral mutation
## rate `u`, which may be zero.
#struct SelectedArchitecture{T,L<:Locus} <: Architecture
#    u    :: T              # mutation rate for neutral loci
#    loci :: Vector{L}      # selected loci
#    xs   :: Vector{Int64}  # location of selected loci
#end
#
## haploid fitness
#function fitness(A::SelectedArchitecture, x::SparseVector{Bool})
#    logw = 0.0 
#    for (i,j) in enumerate(A.xs)
#        # x[1] is either 0 or 1
#        logw += A.loci[i].s * x[j]
#    end
#    exp(logw)
#end
#
#
#"""
#    mutation!(rng, tgt, arch::Architecture)
#"""
#function mutation!(rng::AbstractRNG, tgt, arch::NeutralArchitecture)
#    L = length(tgt)
#    U = arch.u * L
#    n = rand(rng, Poisson(U))
#    n == 0 && return tgt
#    x = sample(rng, 1:L, n, replace=false)
#    @. tgt[x] = !tgt[x]
#    dropzeros!(tgt)  # important!
#    return tgt
#end
#
#"""
#    mutation!(rng, tgt, arch::Architecture)
#"""
#function mutation!(rng::AbstractRNG, tgt, arch::SelectedArchitecture)
#    L = length(tgt)
#    @unpack xs, u, loci = arch
#    # neutral mutations
#    U = arch.u * (L - length(loci))
#    n = rand(rng, Poisson(U))
#    if n != 0
#        x = sample(rng, 1:L, n, replace=false)
#        @. tgt[x] = !tgt[x]
#    end
#    # selected mutations
#    for (i,l) in enumerate(loci)
#        rand(rng) > l.u && continue
#        tgt[xs[i]] = !tgt[xs[i]]
#    end
#    dropzeros!(tgt)  # important!
#    return tgt
#end
#
