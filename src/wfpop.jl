# Implemented for general ploidy levels. See at the bottom of the file for
# specialized functions that need to be implemented.
abstract type Ploid end
struct Haploid <: Ploid end
struct Diploid <: Ploid end
_ploidy(_::Haploid) = 1
_ploidy(_::Diploid) = 2

"""
    WFPopulation

The idea is that at any time, the field `x` has the current haplotypes in
the population.  `_x` serves as a preallocated container in which we store
the offspring while in a generation loop.
"""
@with_kw struct WFPopulation{P<:Ploid,H,T,A<:Architecture,G,R<:RecombinationMap}
    N      :: Int
    arch   :: A 
    gpm    :: G = GPMap()
    recmap :: R 
    ploidy :: P
    nodes  :: Vector{T} = collect(1:_ploidy(ploidy)*N)  # tree sequence nodes
    x      :: Vector{H} = [Bool[] for _=1:_ploidy(ploidy)*N]  # haplotypes
    _x     :: Vector{H} = deepcopy(x) 
end

# Indexing yields an individual's genome
Base.getindex(pop::WFPopulation{Haploid}, i) = pop.x[i]
Base.getindex(pop::WFPopulation{Diploid}, i) = (pop.x[i], pop.x[pop.N + i])
getnodes(pop::WFPopulation{Haploid}, i) = pop.nodes[i]
getnodes(pop::WFPopulation{Diploid}, i) = (pop.nodes[i], pop.nodes[pop.N+i])
getcopy(pop::WFPopulation{Haploid}, i) = copy(pop.x[i])
getcopy(pop::WFPopulation{Diploid}, i) = (copy(pop.x[i]), copy(pop.x[pop.N+i]))

ploidy(pop::WFPopulation) = _ploidy(pop.ploidy)
nhaplotypes(pop::WFPopulation) = pop.N*ploidy(pop)

init_ts(pop::WFPopulation, L; popid=0) = init_ts(nhaplotypes(pop), L, popid=popid)

#eval_fitness(pop::WFPopulation) = map(i->fitness(pop.arch, pop[i]), 1:pop.N)
eval_fitness(pop::WFPopulation) = map(i->exp(phenotype(pop.gpm, pop[i])), 1:pop.N)

function simplify!(pop::WFPopulation, ts::TreeSequence)
    @unpack nodes, N = pop
    sts = simplify(ts, nodes, keep_roots=true)
    nv = length(sts.nodes)
    pop.nodes .= collect(nv-nhaplotypes(pop)+1:nv)
    return pop, sts
end

# no ts recording
function generation!(rng, pop)
    w = eval_fitness(pop)
    idx = sample(rng, 1:pop.N, Weights(w), 2pop.N)
    pop_ = generation!(rng, pop, idx)
end

# with ts recording
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
        generate_offspring!(rng, pop, k, idx[k], idx[N+k], ts, ns)
    end
    _x = mutation!(rng, _x, arch)
    reconstruct(pop, x=_x, _x=x, nodes=collect(ns))
end

function generation!(rng::AbstractRNG, pop::WFPopulation, idx::Vector{Int})
    @unpack N, x, _x, arch, recmap, nodes = pop
    @assert length(idx) == 2N  "Biparental reproduction" 
    for k=1:N  # offspring individual k
        # offspring k has mother and father idx[k] and idx[N+k]
        generate_offspring!(rng, pop, k, idx[k], idx[N+k])
    end
    _x = mutation!(rng, _x, arch)
    reconstruct(pop, x=_x, _x=x)
end

# Functions specialized to ploidy level
function generate_offspring!(rng, pop::WFPopulation{Diploid}, k, p1, p2, args...)
    _generate_offspring!(rng, pop,       k, p1, pop.N+p1, args...)
    _generate_offspring!(rng, pop, pop.N+k, p2, pop.N+p2, args...)
end

function generate_offspring!(rng, pop::WFPopulation{Haploid}, args...)
    _generate_offspring!(rng, pop, args...)
end

function migrate_copy!(src::WFPopulation{Diploid}, dest::WFPopulation{Diploid}, i, k)
    _migrate_copy!(src, dest,       i,       k)
    _migrate_copy!(src, dest, src.N+i, src.N+k)
end

function migrate_copy!(src::WFPopulation{Haploid}, dest::WFPopulation{Haploid}, i, k)
    _migrate_copy!(src, dest,   i,   k)
end

# These are haplotype level functions
# with ts recording
function _generate_offspring!(rng, pop, k, p1, p2, ts::TreeSequence, ns)
    @unpack arch, recmap, nodes = pop 
    (p1, p2) = rand(rng) < 0.5 ? (p1, p2) : (p2, p1)
    bps = rand_breakpoints(rng, recmap)
    recombine!(pop._x[k], bps, pop.x[p1], pop.x[p2], arch.xs) 
    addedges!(ts, nodes[p1], nodes[p2], ns[k], bps)
end

# without ts recording
function _generate_offspring!(rng, pop, k, p1, p2)
    @unpack arch, recmap, nodes = pop 
    (p1, p2) = rand(rng) < 0.5 ? (p1, p2) : (p2, p1)
    bps = rand_breakpoints(rng, recmap)
    recombine!(pop._x[k], bps, pop.x[p1], pop.x[p2], arch.xs) 
end

function _migrate_copy!(src::W, dest::W, i, k) where W<:WFPopulation
    copy!(dest.x[k], src.x[i])
    dest.nodes[k] = src.nodes[i]
    # Note that migration does not change the tree sequence, only which ts
    # nodes are in which population
end

