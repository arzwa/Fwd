abstract type Population end

struct MainlandIsland{I<:Population,M<:Population,T}
    island :: I    # island population
    mainland :: M  # mainland population
    m :: T         # migration rate
end

struct HaploidFixedPopulation{T} <: Population
    genotype :: T
end
nloci(m::HaploidFixedPopulation) = length(m.genotype)

generation!(_, pop::HaploidFixedPopulation) = pop

function rand_migrant!(_::AbstractRNG, x, m::HaploidFixedPopulation) 
    x .= copy(m.genotype)
end

struct HaploidHWLEPopulation{T} <: Population
    p :: T
end
nloci(m::HaploidHWLEPopulation) = length(m.p)

generation!(_, pop::HaploidHWLEPopulation) = pop

function rand_migrant!(rng::AbstractRNG, x, m::HaploidHWLEPopulation) 
    x .= rand(rng, length(x)) .< m.p
end

function rand_migrant(rng::AbstractRNG, m)
    x = Vector{Bool}(undef, nloci(m))
    rand_migrant!(rng, x, m)
end

struct HaploidWFPopulation{T,A,R<:RecombinationMap} <: Population
    x  :: Vector{T}       
    x_ :: Vector{T}
    N  :: Int64
    arch :: A
    recmap :: R
end
nloci(m::HaploidWFPopulation) = length(m.arch)

allele_freqs(pop::HaploidWFPopulation) = sum(pop.x) ./ pop.N 

function rand_migrant!(_::AbstractRNG, x, m::HaploidWFPopulation) 
    x .= copy(rand(m.x))
end

HaploidWFPopulation(x, arch, recmap) = 
    HaploidWFPopulation(x, deepcopy(x), length(x), arch, recmap)
# XXX deepcopy! copy is not safe

function generation!(rng::AbstractRNG,
        metapop::MainlandIsland)
    @unpack mainland, island, m = metapop
    migration!(rng, island, mainland, m)
    island = generation!(rng, island)
    mainland = generation!(rng, mainland)
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
    @inbounds Threads.@threads for k=1:N
        i = idx[k]    # parent 1 index
        j = idx[N+k]  # parent 2 index
        # recombination
        recombination!(x_[k], rng, recmap, x[i], x[j], arch)
        # mutation
        mutations = rand_mutations(rng, arch.mut) 
        mutation!(x_[k], mutations)
    end
    HaploidWFPopulation(x_, x, N, arch, recmap) 
end

