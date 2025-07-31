
@with_kw struct DiploidWFPopulation{H,T,A<:Architecture,R<:RecombinationMap}
    N      :: Int
    x      :: Vector{H}  # haplotypes
    _x     :: Vector{H}  = deepcopy(x) 
    arch   :: A
    recmap :: R 
    nodes  :: Vector{T}  # tree sequence nodes
end

# Indexing yields an individual's genome
Base.getindex(d::DiploidWFPopulation, i) = (d.x[i], d.x[d.N + i])

function eval_fitness(pop::DiploidWFPopulation)
    @unpack N, x, arch = pop
    w = map(i->fitness(arch, x[i], x[N+i]), 1:N)
end

function simplify!(pop::DiploidWFPopulation, ts::TreeSequence)
    @unpack nodes, N = pop
    sts = simplify(ts, nodes, keep_roots=true)
    nv = length(sts.nodes)
    pop.nodes .= collect(nv-2N+1:nv)
    return pop, sts
end

"""
`idx` is a vector of length 2N with numbers ∈ [1..N]  where entry `k` and
entry `N+k` correspond to the indices of the mother and father of offspring
individual `k`. This means that the maternal haplotypes are at `idx[k]` and
`N+idx[k]`, and the paternal haplotypes at `idx[N+k]` and `N+idx[N+k]`.
"""
function generation!(
        rng::AbstractRNG,
        pop::DiploidWFPopulation, 
        idx::Vector{Int},
        ts::TreeSequence,
        popid=0)   # XXX don't like the `popid`
    @unpack N, x, _x, arch, recmap, nodes = pop
    @assert length(idx) == 2N
    # new nodes to ts
    exnode = ts.nodes[nodes[1]]
    ns = addnodes!(ts, 2N, Node(time(exnode)+1, popid))
    for k=1:N  # offspring individual k
        _generate_offspring!(rng, pop, ts, idx[k], ns[k], k)
        _generate_offspring!(rng, pop, ts, idx[N+k], ns[N+k], N+k)
    end
    DiploidWFPopulation(N, _x, x, arch, recmap, collect(ns))
end

function _generate_offspring!(rng, pop, ts, p, o, i)
    @unpack N, arch, recmap, nodes = pop 
    (p1, p2) = rand(rng) < 0.5 ? (p, N+p) : (N+p, p)
    bps, edges = recombine(rng, nodes[p1], nodes[p2], o, recmap)
    recombine!(pop._x[i], bps, pop.x[p1], pop.x[p2], arch.xs) 
    addedges!(ts, edges)
    mutations = rand_mutations(rng, arch.mut) 
    mutation!(pop._x[i], mutations)
end


"""
    TwoPopOneWay

Migration is from A to B forward in time. Migration replaces an expected
proportion `m` of the B population by A individuals, without affecting A.
"""
struct TwoPopOneWay{P,T}
    m    :: T
    popA :: P
    popB :: P
end

active_nodes(m::TwoPopOneWay) = [m.popA.nodes ; m.popB.nodes]

"""
    migration!(rng, metapop)

!!! This migration function implements migration-as-copying, i.e. a
proportion `m` of the B population is replaced by indviduals from A, but A
is unaffected.  (biologically this could correspond to sending out asexual
propagules -- although this is more relevant for haplodiplontic cryptogams
etc. than for diploids).
"""
function migration!(rng, metapop::TwoPopOneWay{P}) where P<:DiploidWFPopulation
    @unpack popA, popB, m = metapop
    NB = popB.N
    NA = popA.N
    nmig = min(NB, rand(rng, Poisson(m*NB)))
    idx = sample(rng, 1:NA, nmig) 
    for k=1:nmig
        i = idx[k]  # i is the index of the migrant individual in A
        copy!(popB.x[k]   , popA.x[i])
        copy!(popB.x[NB+k], popA.x[NA+i])
        popB.nodes[k]    = popA.nodes[i]
        popB.nodes[NB+k] = popA.nodes[NA+i]
    end
end

# Other possibilities for migration:
# 1. as in SLiM, a proportion m of the offspring is replaced by offspring
# from a pair of individuals from the source population
# 2. ...

function generation!(rng, metapop::TwoPopOneWay, ts::TreeSequence)
    migration!(rng, metapop)
    @unpack popA, popB = metapop
    w = eval_fitness(popA)
    idx = sample(rng, 1:popA.N, Weights(w), 2popA.N)
    popA_ = generation!(rng, popA, idx, ts, 1)
    w = eval_fitness(popB)
    idx = sample(rng, 1:popB.N, Weights(w), 2popB.N)
    popB_ = generation!(rng, popB, idx, ts, 2)
    reconstruct(metapop, popA=popA_, popB=popB_)
end

function simplify!(pop::TwoPopOneWay, ts::TreeSequence)
    @unpack popA, popB = pop
    ns = active_nodes(pop)
    sts = simplify(ts, ns, keep_roots=true)
    #T = time(sts.nodes[end])
    # XXX the indices are not what I expected?
    #popA.nodes .= findall(x->time(x) == T && population(x) == 1, sts.nodes)
    #popB.nodes .= findall(x->time(x) == T && population(x) == 2, sts.nodes)
    nv = length(sts.nodes)
    NA = length(popA.nodes)
    NB = length(popB.nodes)
    na = (nv-NA-NB+1):(nv-NB)
    nb = (nv-NB+1):nv
    popA.nodes .= collect(na)
    popB.nodes .= collect(nb)
    return pop, sts
end

function init_ts(pop::TwoPopOneWay{P}, L::V) where {P<:DiploidWFPopulation,V}
    @unpack popA, popB = pop
    N = popA.N + popB.N
    nodesA = [Node(0, 1) for _=1:2popA.N]
    nodesB = [Node(0, 2) for _=1:2popB.N]
    edges = Edge{Int,V}[]
    children = [Int[] for _=1:2N]
    TreeSequence([nodesA; nodesB], edges, children, L, true)
end

