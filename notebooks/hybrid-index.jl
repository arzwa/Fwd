# Hybrid index, or rather proportion of migrant ancestry simulations.
# This is not that straightforward.
# One thing one could do is, at time `t`, sample some individual and look
# for the first generation in the past where it had a migrant ancestor in
# the pedigree, then look down from that migrant ancestor how many
# offspring it had one, two, ... generations later and how much genetic
# material it left for each of these. We would need an *unsimplified* ts
# for this.
# Actually, one need not start from a contemporary individual. We could
# just simulate a population, and before simplifying the tree sequence at
# some regular interval scan the `ts` starting from the previous
# simplification time point and look for migrants and trace their ancestry. 

# I'll start from the same simulation as in `divsel-polygenic-2.jl`
# (2026-03-23)

function sim(C)
    rng = Random.seed!(135)
    Ls  = 0.25
    L   = 50
    #Ls  = 0.1
    #L   = 20
    s̄   = Ls/L
    dfe = Exponential(s̄)
    ss  = rand(rng, dfe, L)
    ss .*= s̄/mean(ss) 
    α   = 2.0
    zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
    ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
    u   = s̄/200
    NA  = 1
    Ns  = 5. 
    NB  = ceil(Int64, Ns/s̄)
    ms  = 0.3
    m   = ms*s̄
    xs  = ys .* C
    BP  = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u), α=0.2)
    R   = LinearMap(C)
    AA  = Architecture([BiAllelic(0.0)   for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u)     for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
    MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
    nA  = collect(1:NA)
    nB  = collect(1:NB) .+ NA
    xA  = [ ones(Int, L) for _=1:NA]
    xB  = [zeros(Int, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    ngen = 50NB
    mpop, qs = Fwd.simulate!(mpop, ngen, pop->mean(pop.popB.x), every=100)
    # we can start with a new tree sequence
    pop = deepcopy(mpop)
    nA  = collect(1:NA)
    nB  = collect(1:NB) .+ NA
    pop.popA.nodes .= nA
    pop.popB.nodes .= nB
    ts = init_ts(pop)
    ngen2 = 5000
    pop, ts, qs_ = Fwd.simulate!(pop, ts, ngen2, 
        pop->mean(pop.popB.x), simplify=Inf)
    ts, BP, qs
end

maplens = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
ress = map(sim, maplens)

function rho(BP, kmax=floor(Int, log2(0.5/BP.model.m)))
    # heuristic for kmax: at most 50% non-residents
    @unpack m, s, xs = BP.model
    ws = map(0:kmax-1) do k
        exp(-sum([s[i]*BP.Ep[i] for i=1:length(xs)])/2^k)
    end 
    fs = m*cumprod(2 .* ws)
    ρ = sum([1/2^k * fs[k] for k=1:kmax])
    ρ, fs, ws
end

# Now generalize. This is rather complicated!
function migrant_ancestries(ts, kmax, NB, C)
    tmax = ts[end].time
    # mainland individuals
    Mnodes = filter(i->ts[i].pop == 1 && ts[i].time != tmax, 1:length(ts))
    res = trace_ancestry_onegen(
        Mnodes, [[(0.,C)] for i=1:length(Mnodes)], ts)
    results = [res]
    Cs = zeros(tmax, kmax)
    k = 1
    while k <= kmax
        # take those that did inherit migrant ancestry
        res = filter(x->!isempty(last(x)), res)
        for (_, t, off) in res
            Cs[t+1,k] += sum(last.(off)) / (NB*C)
        end 
        X = vcat(map(x->x[3], res)...)
        k += 1
        nodes = first.(X)
        blocks = getindex.(X,2)
        res = trace_ancestry_onegen(nodes, blocks, ts)
        push!(results, res)
    end
    return Cs, results
end

"""
    trace_ancestry_onegen
For a set of parent nodes, with associated blocks of genome, obtain
all descendants and the intersection of what they inherited from the parent 
and the blocks associated with that parent.
"""
function trace_ancestry_onegen(nodes, blocks, ts)
    map(zip(nodes, blocks)) do (p, pblocks)
        pedges = ts.adjlist[p]
        children = unique(map(e->ts.edges[e].child, pedges))
        BCk = filter(c->ts[c].pop == 2, children)
        off = map(BCk) do c
            cblocks = [(e.left, e.rght) for 
                e in ts.edges[pedges] if e.child == c]
            oblocks = intersect_blocks(pblocks, cblocks)
            spans = sum(map(x->x[2]-x[1], oblocks))
            c, oblocks, spans
        end
        p, ts[p].time, off
    end
end

function intersect_blocks(pblocks::Vector{T}, cblocks::Vector{T}) where T
    blocks = T[]
    for (x0, x1) in pblocks
        for (y0, y1) in cblocks
            y0 > x1 && continue
            y1 < x0 && continue
            z0 = y0 > x0 ? y0 : x0
            z1 = x1 < y1 ? x1 : y1
            push!(blocks, (z0, z1))
        end
    end
    return blocks
end

function avr(xs, t)
    RM = Fwd.rec_matrix(xs)
    RT = (1 .- RM) .^ t
    L = length(xs)
    r̄ = 0.
    for i=2:L
        for j=1:i-1
            r̄ += RT[i,j]
        end
    end
    r̄ *= 2/(L*(L-1))
end

function rpreds(BP, kmax)
    @unpack m, s, xs = BP.model
    _, _, ws = rho(BP, kmax)
    z = m*ws[1]
    S = sum(BP.Ep .* s)
    zs = [z]
    for k=1:kmax-1
        #z *= (1 - 0.5S*avr(xs, k))
        #z *= exp(-0.5S*avr(xs, k))
        z *= exp(-S*avr(xs, k))
        push!(zs, z)
    end
    zs
end

kmax = 10 # XXX blows up for k large (e.g. kmax=20)
anc = map(zip(maplens, ress)) do (C, (ts, BP, _))
    Cs, res = migrant_ancestries(ts, kmax, BP.model.Ne, C)
    cs = mean(Cs, dims=1) |> vec
    fs, ws = Barriers.bc_fitnesses(BP, kmax)
    unlinked = [fs[k]/2^k for k=1:kmax]
    linked = Barriers.bc_ancestries(BP, kmax)
    C, cs, linked, unlinked
end 

PP = map(enumerate(anc)) do (i,(C, cs, linked, unlinked))
    plot(unlinked, label="unlinked")
    plot!(linked, label="linkage")
    scatter!(cs, marker=true, ms=2, label="", legend=:topright,
        title="$C M", ylabel=i ∈ [1,4] ? "prop. of population\nw/ migrant ancestry" : "", 
        color=:black, xlabel=i > 3 ? "generation" : "")
    xticks!(1:2:kmax, [["F1"];  ["BC$k" for k=2:2:kmax-1]])
end
plot(PP..., layout=(2,3), size=(700,400))


# number of outgoing edges
nout = map(y->map(x->length(x[3]), y), res)
# For short map lengths (e.g. looking at a single locus), we should have
# mean(nout[1]) ≈ m*ws[1]*NB.
# For long map lengths, we should have mean(nout[1]) ≈ 2m*ws[1]*NB
# the rest should be in between. For short maps, the variance should be
# twice the mean.
# We do see this if we do a neutral simulation with C = 1e-14 for instance.

C = 1.0
exts, exbp, qs = sim(C)

kmax = 15
Cs, res = migrant_ancestries(exts, kmax, exbp.model.Ne, C)
cs = mean(Cs, dims=1) |> vec
fs, ws = Barriers.bc_fitnesses(exbp, kmax)
unlinked = [fs[k]/2^k for k=1:kmax]
linked = Barriers.bc_ancestries(exbp, kmax)

plot(unlinked, label="unlinked")
plot!(linked, label="linkage")
scatter!(cs, marker=true, ms=2, label="", legend=:topright, color=:black)
xticks!(1:2:kmax, [["F1"];  ["BC$k" for k=2:2:kmax-1]])
