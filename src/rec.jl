abstract type RecombinationMap end

struct LinearMap{T} <: RecombinationMap
    maplength :: T  # maplength in Morgans, i.e. expected # of crossovers
end
Base.length(m::LinearMap) = m.maplength
Base.zero(m::LinearMap{T}) where T = zero(T)
function nextbreakpoint(rng, recmap::LinearMap, bp)
    @unpack maplength = recmap
    maplength == 0. && return maplength
    bp = bp + randexp(rng)
    return min(bp, maplength)
end

@with_kw struct LinearPhysMap{T} <: RecombinationMap
    G :: Int
    C :: T
    rbp :: T = C/(G-1)   # M/bp
end
Base.length(m::LinearPhysMap) = m.G
Base.zero(m::LinearPhysMap) = 0
recrate(m::LinearPhysMap, x, y) = recrate(m.rbp * abs(x - y))
function nextbreakpoint(rng, recmap::LinearPhysMap, bp)
    @unpack rbp, G = recmap
    (isnan(rbp) || iszero(rbp)) && return G
    bp = bp + ceil(Int, rand(rng, Exponential(1/rbp)))
    return min(bp, G)
end

struct Unlinked <: RecombinationMap end
Base.length(m::Unlinked) = 1
Base.zero(m::Unlinked) = 0
nextbreakpoint(_, recmap::Unlinked, bp) = 1

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
M = Chromosomes([LinearPhysMap(G=1, C=0.0) for _=1:L])
```
The latter however admits mixing linked with unlinked stuff.
"""
struct Chromosomes{M<:RecombinationMap} <: RecombinationMap 
    maps :: Vector{M}
end    
Base.length(c::Chromosomes) = length(c.maps)

# Haldane's mapping function
# distance -> recombination rate
recrate(d) = 0.5*(1-exp(-2d))     

# recombination rate -> distance
distance(r) = -0.5*log(1 - 2r)  

# recombination rate matrix
rec_matrix(x) = [recrate(abs(x[i] - x[j])) for i=1:length(x), j=1:length(x)]

# Assumes each `recmap` implements `length` and `nextbreakpoint`
# length gives the map lengths relative to which `xs` are coordinates.
function recombine!(rng, tgt, src1, src2, recmap::Chromosomes, xs, args...)
    χ  = true  # indicator which src to take
    i  = 1     # locus
    C  = 0.0   # map position start
    for (k,chrom) in enumerate(recmap.maps)
        i, C = recombine!(rng, 
            tgt, src1, src2, chrom, xs, args...; 
            i=i, C=C, χ=χ)
        χ = rand(rng) < 0.5   # recombination between chromosomes
    end
end

# without ts recording
function recombine!(rng, tgt, src1, src2, recmap, xs; 
        i=1, C=zero(recmap), χ=true)
    bp = zero(recmap)
    C′ = C + length(recmap)
    while C + bp < C′
        bp = nextbreakpoint(rng, recmap, bp)
        bpk = C + bp  # absolute breakpoint coordinate 
        # if bp is at `xs[i]`, than xs[i] is the last locus before the bp
        while i <= length(xs) && xs[i] <= bpk && xs[i] <= C′
            tgt[i] = χ ? src1[i] : src2[i]
            i += 1
        end  
        χ = !χ  # switch parent
    end
    return i, C′
end

# with ts recording
function recombine!(rng, tgt, src1, src2, recmap, xs, ts, nodes; 
        i=1, C=zero(recmap), χ=true)
    (p1, p2, c) = nodes
    bp = zero(recmap)
    x0 = C
    C′ = C + length(recmap)
    while x0 < C′
        bp = nextbreakpoint(rng, recmap, bp)
        x1 = C + bp
        e  = χ ? TS.Edge(p1, c, x0, x1) : TS.Edge(p2, c, x0, x1)
        TS.addedge!(ts, e) 
        # if bp is at `xs[i]`, than xs[i] is the last locus before the bp
        while i <= length(xs) && xs[i] <= x1 && xs[i] <= C′
            tgt[i] = χ ? src1[i] : src2[i]
            i += 1
        end  
        x0 = x1 
        χ  = !χ  # switch parent
    end
    return i, C
end

