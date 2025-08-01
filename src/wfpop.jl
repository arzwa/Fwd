# Implemented for general ploidy levels. See at the bottom of the file for
# specialized functions that need to be implemented.
abstract type Ploid end
struct Haploid <: Ploid end
struct Diploid <: Ploid end
ploidy(_::Haploid) = 1
ploidy(_::Diploid) = 2

"""
    WFPopulation

The idea is that at any time, the field `x` has the current haplotypes in
the population.  `_x` serves as a preallocated container in which we store
the offspring while in a generation loop.
"""
@with_kw struct WFPopulation{P<:Ploid,H,T,A<:Architecture,R<:RecombinationMap}
    N      :: Int
    x      :: Vector{H}  # haplotypes
    _x     :: Vector{H}  = deepcopy(x) 
    arch   :: A
    recmap :: R 
    nodes  :: Vector{T}  # tree sequence nodes
    ploidy :: P
end

# Indexing yields an individual's genome
Base.getindex(pop::WFPopulation{Haploid}, i) = pop.x[i]
Base.getindex(pop::WFPopulation{Diploid}, i) = (pop.x[i], pop.x[pop.N + i])

ploidy(pop::WFPopulation) = ploidy(pop.ploidy)
nhaplotypes(pop::WFPopulation) = pop.N*ploidy(pop)

init_ts(pop::WFPopulation, L; popid=0) = init_ts(nhaplotypes(pop), L, popid=popid)

eval_fitness(pop::WFPopulation) = map(i->fitness(pop.arch, pop[i]), 1:pop.N)

function simplify!(pop::WFPopulation, ts::TreeSequence)
    @unpack nodes, N = pop
    sts = simplify(ts, nodes, keep_roots=true)
    nv = length(sts.nodes)
    pop.nodes .= collect(nv-nhaplotypes(pop)+1:nv)
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
        pop::WFPopulation, 
        idx::Vector{Int},
        ts::TreeSequence,
        popid=0)   # XXX don't like the `popid`
    @unpack N, x, _x, arch, recmap, nodes = pop
    @assert length(idx) == 2N  "Biparental reproduction" 
    # new nodes to ts
    exnode = ts.nodes[nodes[1]]
    ns = addnodes!(ts, length(nodes), Node(time(exnode)+1, popid))
    for k=1:N  # offspring individual k
        # offspring k has mother and father idx[k] and idx[N+k]
        generate_offspring!(rng, pop, ts, ns, k, idx[k], idx[N+k])
    end
    reconstruct(pop, x=_x, _x=x, nodes=collect(ns))
end

# These are hplotype level functions
function _generate_offspring!(rng, pop, ts, ns, k, p1, p2)
    @unpack arch, recmap, nodes = pop 
    (p1, p2) = rand(rng) < 0.5 ? (p1, p2) : (p2, p1)
    bps, edges = recombine(rng, nodes[p1], nodes[p2], ns[k], recmap)
    recombine!(pop._x[k], bps, pop.x[p1], pop.x[p2], arch.xs) 
    addedges!(ts, edges)
    mutations = rand_mutations(rng, arch.mut) 
    mutation!(pop._x[k], mutations)
end

function _migrate!(src::W, dest::W, i, k) where W<:WFPopulation
    copy!(dest.x[k] , src.x[i])
    dest.nodes[k] = src.nodes[i]
end

# Functions specialized to ploidy level
function generate_offspring!(rng, pop::WFPopulation{Diploid}, ts, ns, k, p1, p2)
    _generate_offspring!(rng, pop, ts, ns,       k, p1, pop.N+p1)
    _generate_offspring!(rng, pop, ts, ns, pop.N+k, p2, pop.N+p2)
end

function generate_offspring!(rng, pop::WFPopulation{Haploid}, args...)
    _generate_offspring!(rng, pop, args...)
end

function migrate!(src::WFPopulation{Diploid}, dest::WFPopulation{Diploid}, i, k)
    _migrate!(src, dest,   i,   k)
    _migrate!(src, dest, N+i, N+k)
end

function migrate!(src::WFPopulation{Haploid}, dest::WFPopulation{Haploid}, i, k)
    _migrate!(src, dest,   i,   k)
end


