
# metapopulation with migration as copying, that is, m₁₂ is the expected
# proportion of migrants from 1 into 2 (forward in time), meaning that an
# expected m₁₂ individuals from population 2 are replaced by clones copied
# from 1.
struct MetaPop{Pop,T} <: AbstractPop
    P :: Vector{Pop}
    M :: Matrix{T} 
end
Base.length(m::MetaPop) = length(m.P)
Base.getindex(m::MetaPop, i) = m.P[i]
active_nodes(m::MetaPop) = mapreduce(p->p.nodes, vcat, m.P)

struct Migrant{T,V}
    pop :: Int
    idx :: Int
    nodes :: T
    haplotypes :: V
end

_migrants(metapop::MetaPop{<:WFPopulation{Haploid,X}}) where X = Migrant{Int,X}[]
_migrants(metapop::MetaPop{<:WFPopulation{Diploid,X}}) where X = Migrant{Tuple{Int,Int},Tuple{X,X}}[]

# With TwoPopOneWay we didn't have to worry much about how to implement
# migration. Here we have to (we cannot naively do mig. 1->2 and then 2->1,
# since that would allow back-and-forth migration in one gen etc.)
function migration!(rng, metapop)
    @unpack P, M = metapop
    # We have to first collect the migrants ...
    migrants = _migrants(metapop)
    for i=1:length(P)
        # migration into i
        Ni = P[i].N
        k0 = 1
        for j=1:length(P)
            i == j && continue
            M[j,i] == 0. && continue
            Nj = P[j].N
            nmig = min(Ni, rand(rng, Poisson(M[j,i]*Ni)))
            idx = sample(rng, 1:Nj, nmig) 
            for (l,k) = enumerate(k0:k0+nmig-1)
                # k is the index the migrant will get in i
                # idx[l] is the index of the migrant source individual in j 
                mig = Migrant(i, k, getnodes(P[j], idx[l]), getcopy(P[j], idx[l]))
                push!(migrants, mig) 
            end
            k0 += nmig
        end
    end
    # ... and then execute migration.
    for mig in migrants
        _migrate!(P[mig.pop], mig)
    end
    return metapop
end

function _migrate!(pop::WFPopulation{Haploid}, migrant)
    @unpack idx, nodes, haplotypes = migrant
    pop.nodes[idx] = nodes
    pop.x[idx] = haplotypes
end

function _migrate!(pop::WFPopulation{Diploid}, migrant)
    @unpack idx, nodes, haplotypes = migrant
    pop.nodes[idx] = nodes[1]
    pop.nodes[idx+pop.N] = nodes[2]
    pop.x[idx] = haplotypes[1]
    pop.x[idx+pop.N] = haplotypes[2]
end

function generation!(rng, metapop::MetaPop, ts::TreeSequence)
    metapop = migration!(rng, metapop)
    _P = map(1:length(metapop)) do i
        generation!(rng, metapop[i], ts, popid=i)
    end
    reconstruct(metapop, P=_P)
end

function generation!(rng, metapop::MetaPop)
    metapop = migration!(rng, metapop)
    _P = map(1:length(metapop)) do i
        generation!(rng, metapop[i])
    end
    reconstruct(metapop, P=_P)
end

function init_ts(pop::MetaPop)
    tss = map(i->init_ts(pop[i], popid=i), 1:length(pop))
    ns = mapreduce(ts->ts.nodes, vcat, tss)
    es = tss[1].edges  # empty anyhow
    cs = mapreduce(ts->ts.adjlist, vcat, tss)
    TreeSequence(ns, es, cs, tss[1].L, true)
end

function simplify!(pop::MetaPop, ts::TreeSequence)
    ns = active_nodes(pop)
    sts = TS.simplify(ts, ns, keep_roots=true)
    nv = length(sts.nodes)
    for i=length(pop):-1:1
        Ni = length(pop[i].nodes)
        ni = nv-Ni+1:nv
        pop[i].nodes .= collect(ni)
        nv -= Ni
    end
    return pop, sts
end

# split a population in two: sample NA individuals with replacement from
# pop to constitute popA, NB to constitute popB. Update the ts accordingly.
# This uses the metapop migration implementation (we migrate from the
# ancestral in the new subpopulations.
function splitpop(pop, ts, Ns; popids=(1,2))
    pops = map(zip(Ns, popids)) do (N, i)
        n = _ploidy(pop.ploidy)*N
        pop_ = reconstruct(pop, 
            x  = similar(pop.x, n), 
            _x = similar(pop.x, n),
            nodes = similar(pop.nodes, n),
            N = N)
        idx = sample(1:pop.N, N, replace=true) 
        for k in 1:N
            mig = Migrant(i, k, 
                getnodes(pop, idx[k]), 
                getcopy(pop, idx[k]))
            _migrate!(pop_, mig)
        end
        pop_ = reconstruct(pop_, _x=copy.(pop_.x))
        pop_
    end
end
