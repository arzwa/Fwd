
# metapopulation with migration as copying, that is, m₁₂ is the expected
# proportion of migrants from 1 into 2 (forward in time), meaning that an
# expected m₁₂ individuals from population 2 are replaced by clones copied
# from 1.
struct MetaPop{Pop,T}
    P :: Vector{Pop}
    M :: Matrix{T} 
end
Base.length(m::MetaPop) = length(m.P)
Base.getindex(m::MetaPop, i) = m.P[i]
active_nodes(m::MetaPop) = mapreduce(p->p.nodes, vcat, m.P)

# With TwoPopOneWay we didn't have to worry much about how to implement
# migration. Here we have to (we cannot naively do mig. 1->2 and then 2->1,
# since that would allow back-and-forth migration in one gen etc.)
# We work x -> _x, and after all migration has been determined, reconstruct
# the model x=_x, _x=x.
function migration!(rng, metapop::MetaPop)
    @unpack P, M = metapop
    bookkeeping = NTuple{3,Int64}[]
    for i = 1:length(P)
        # migration into i
        Ni = P[i].N
        for j=1:length(P)
            i == j && continue
            M[j,i] == 0. && continue
            Nj = P[j].N
            nmig = min(Ni, rand(rng, Poisson(M[j,i]*Ni)))
            idx = sample(rng, 1:Nj, nmig) 
            for k=1:nmig
                push!(bookkeeping, (i, k, P[j].nodes[idx[k]]))
            end
        end
    end
    for (i, k, nk) in bookkeeping
        P[i].nodes[k] = nk
    end
    # ... how to do this best?
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

function init_ts(pop::MetaPop, L)
    tss = map(i->init_ts(pop[i], L, popid=i), 1:length(pop))
    ns = mapreduce(ts->ts.nodes, vcat, tss)
    es = tss[1].edges  # empty anyhow
    cs = mapreduce(ts->ts.children, vcat, tss)
    TreeSequence(ns, es, cs, L, true)
end

function simplify!(pop::MetaPop, ts::TreeSequence)
    ns = active_nodes(pop)
    sts = simplify(ts, ns, keep_roots=true)
    nv = length(sts.nodes)
    for i=length(pop):-1:1
        Ni = length(pop[i].nodes)
        ni = nv-Ni+1:nv
        pop[i].nodes .= collect(ni)
        nv -= Ni
    end
    return pop, sts
end
