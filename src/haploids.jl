abstract type Population end

struct MainlandIsland{I<:Population,M<:Population,T}
    island :: I    # island population
    mainland :: M  # mainland population
    m :: T         # migration rate
end

struct HaploidFixedPopulation{T} <: Population
    genotype :: T
end

function rand_migrant!(_::AbstractRNG, x, m::HaploidFixedPopulation) 
    x .= copy(m.genotype)
end

struct HaploidHWLEPopulation{T} <: Population
    p :: T
end

function rand_migrant!(rng::AbstractRNG, x, m::HaploidHWLEPopulation) 
    x .= rand(rng, length(x)) .< m.p
end

function rand_migrant(rng::AbstractRNG, m::HaploidHWLEPopulation)
    x = Vector{Bool}(undef, length(m.p))
    rand_migrant!(rng, x, m)
end

struct HaploidWFPopulation{T,A,R} <: Population
    x  :: Vector{T}       
    x_ :: Vector{T}
    N  :: Int64
    arch :: A
    recmap :: R
end

HaploidWFPopulation(x, arch, recmap) = 
    HaploidWFPopulation(x, deepcopy(x), length(x), arch, recmap)
# XXX deepcopy! copy is not safe

function generation!(rng::AbstractRNG,
        metapop::MainlandIsland)
    @unpack mainland, island, m = metapop
    migration!(rng, island, mainland, m)
    island = generation!(rng, island)
    return MainlandIsland(island, mainland, m)
end

function migration!(rng::AbstractRNG,
        island::HaploidWFPopulation, 
        mainland, m)
    @unpack N = island
    nmig = min(N, rand(rng, Poisson(m*N)))
    for i=1:nmig
        rand_migrant!(rng, island.x[i], mainland)
    end
end

function generation!(rng::AbstractRNG, 
        pop::HaploidWFPopulation)
    @unpack N, x, x_, arch, recmap = pop
    w   = map(xi->fitness(arch, xi), x)
    idx = sample(rng, 1:N, Weights(w), 2N)
    @inbounds for k=1:N
        i = idx[k]    # parent 1 index
        j = idx[N+k]  # parent 2 index
        # recombination
        breakpoints = rand_breakpoints(rng, recmap)
        recombine!(x_[k], breakpoints, x[i], x[j], arch.xs)
        # mutation
        mutations = rand_mutations(rng, arch.mut) 
        mutation!(x_[k], mutations)
    end
    HaploidWFPopulation(x_, x, N, arch, recmap) 
end

