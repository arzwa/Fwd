
@with_kw struct DiploidWFPopulation{H,T,A<:Architecture,R<:RecombinationMap}
    N      :: Int
    x      :: Vector{H}  # haplotypes
    _x     :: Vector{H}  = deepcopy(x) 
    arch   :: A
    recmap :: R 
    nodes  :: Vector{T}  # tree sequence nodes
end

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
        m = idx[k]
        f = idx[N+k]
        om = ns[k]
        of = ns[N+k]
        # maternal haplotype
        (m1, m2) = rand(rng) < 0.5 ? (m, N+m) : (N+m, m)
        bps, edges = recombine(rng, nodes[m1], nodes[m2], om, recmap)
        recombine!(_x[k], bps, x[m1], x[m2], arch.xs) 
        addedges!(ts, edges)
        if length(bps) > 1
            nk, n1, n2 = om, nodes[m1], nodes[m2]
           # @info "" (nk, n1, n2) (k, _x[k]) (m1, x[m1]) (m2, x[m2]) bps edges
        end
        mutations = rand_mutations(rng, arch.mut) 
        mutation!(_x[k], mutations)
        # paternal haplotype
        (f1, f2) = rand(rng) < 0.5 ? (f, N+f) : (N+f, f)
        bps, edges = recombine(rng, nodes[f1], nodes[f2], of, recmap)
        recombine!(_x[N+k], bps, x[f1], x[f2], arch.xs) 
        addedges!(ts, edges)
        if length(bps) > 1
            nk, n1, n2 = of, nodes[f1], nodes[f2]
          #  @info "" (nk, n1, n2) (k, _x[k]) (f1, x[f1]) (f2, x[f2]) bps edges
        end
        mutations = rand_mutations(rng, arch.mut) 
        mutation!(_x[N+k], mutations)
        # XXX better have mutations at population level, since arch is
        # population level anyhow?
    end
    DiploidWFPopulation(N, _x, x, arch, recmap, collect(ns))
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
        #copy!(popB.x[k]   , copy(popA.x[i]))
        #copy!(popB.x[NB+k], copy(popA.x[NA+i]))
        #XXX paranoid
        popB.x[k]    = deepcopy(popA.x[i])
        popB.x[NB+k] = deepcopy(popA.x[NA+i])
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
    sts = simplify(ts, ns)
    T = time(sts.nodes[end])
    # XXX the indices are not what I expected?
    popB.nodes .= findall(x->time(x) == T && population(x) == 2, sts.nodes)
    popA.nodes .= findall(x->time(x) == T && population(x) == 1, sts.nodes)
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

