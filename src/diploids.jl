
struct DiploidWFPopulation{H,A<:Architecture,R<:RecombinationMap}
    N  :: Int
    arch :: A
    recmap :: R 
    x  :: Vector{H}  # haplotypes
    x_ :: Vector{H}  
end

function generation!(rng::AbstractRNG, 
        pop::DiploidWFPopulation, 
        tsrec::TreeSequenceRecorder)
    @unpack N, x, x_, arch, recmap = pop
    @unpack active, ts = tsrec
    # fitness calculation -- XXX could be taken out to have a function that
    # only deals with reproduction...
    w = map(i->fitness(arch, x[i], x[N+i]), 1:N)
    # sample 2N parents
    idx = sample(rng, 1:N, Weights(w), 2N)
    # new nodes to ts
    ns = length(ts.nodes)+1:length(ts.nodes)+2N 
    for k=1:N  # offspring individual k
        for (i, p, o) in [(k, idx[k], ns[k]), (N+k, idx[N+k], ns[N+k])]
            # `i` is the index of the offspring haplotype we're creating
            # `p` is the index of the parent individual whose haplotypes we recombine
            # `o` is the node ID in the tree sequence of the offspring haplotype
            bps, edges = recombine(rng, active[p], active[N+p], o, recmap)
            for e in edges
                addedge!(ts.edges, ts.children, e)
            end
            recombine!(x_[i], bps, x[p], x[N+p], arch.xs) 
            mutations = rand_mutations(rng, arch.mut) 
            mutation!(x_[i], mutations)
            # XXX better have mutations at population level, since arch is
            # population level anyhow?
        end
    end
    addnodes!(ts, 2N, ts.nodes[active[1]]+1)
    tsrec.active .= ns
    DiploidWFPopulation(N, arch, recmap, x_, x)
end

function recombine(rng, p1, p2, c, recmap::RecombinationMap)
    bps = rand_breakpoints(rng, recmap)
    edges = map(enumerate(bps)) do (i,x1)
        x0 = i == 1 ? zero(x1) : bps[i-1]
        isodd(i) ? Edge(p1, c, x0, x1) : Edge(p2, c, x0, x1)
    end
    bps, edges
end

function TwoPopUni{P,T}
    m :: T
    popA :: P
    popB :: P
end

function migration!(rng, metapop::TwoPopUni{P}) where P<:DiploidWFPopulation
    @unpack N = metapop.B
    nmig = min(N, rand(rng, Poisson(m*N)))
    for i=1:nmig
        rand_migrant!(rng, metapop.popB, metapop.popA, i)
    end
end

function rand_migrant!(rng, popB, popA, i)
    k = rand(rng, 1:popA.N)
    copy!(popB.x[i], popA.x[k])
    copy!(popB.x[i+N], popA.x[k+N])
end
