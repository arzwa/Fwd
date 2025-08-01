"""
    DiploidWFPopulation

The idea is that at any time, the field `x` has the current haplotypes in
the population.  `_x` serves as a preallocated container in which we store
the offspring while in a generation loop.
"""
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

init_ts(pop::DiploidWFPopulation, L; popid=0) = init_ts(2pop.N, L, popid=popid)

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

function generation!(rng, pop, ts; popid=0)
    w = eval_fitness(pop)
    idx = sample(rng, 1:pop.N, Weights(w), 2pop.N)
    pop_ = generation!(rng, pop, idx, ts, popid)
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

function migrate!(src::DiploidWFPopulation, dest::DiploidWFPopulation, i, k)
    copy!(dest.x[k] , src.x[i])
    dest.nodes[k] = src.nodes[i]
    copy!(dest.x[dest.N+k], src.x[src.N+i])
    dest.nodes[dest.N+k] = src.nodes[src.N+i]
end

