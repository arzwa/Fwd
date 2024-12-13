
abstract type Population end

struct DiploidWFPopulation{A<:Architecture,R<:RecombinationMap,T} <: Population
    x  :: Vector{T}
    x_ :: Vector{T}
    N  :: Int64
    arch :: A
    recmap :: R
end
    
function DiploidWFPopulation(x, arch, recmap)
    DiploidWFPopulation(x, similar(x), length(x) ÷ 2, arch, recmap)
end

allelefreqs(pop::DiploidWFPopulation) = sum(pop.x) ./ 2pop.N

function Base.show(io::IO, pop::DiploidWFPopulation)
    @unpack x, N, arch, recmap = pop
    write(io, "DiploidWFPopulation\n  $arch\n  $recmap\n  " *
        "L = $(length(x[1])/Mb)Mb\n  N = $N")
end

function generation!(rng::AbstractRNG, 
        pop::DiploidWFPopulation{<:NeutralArchitecture})
    @unpack N, x, x_, arch, recmap = pop
    idx = rand(rng, 1:N, 2N)
    for k=1:N
        i = idx[k]    # parent 1 index
        j = idx[N+k]  # parent 2 index
        # first haplotype of indiv. k is the result of recombination of two
        # haplotypes of individual i (stored at indices i and N+i)
        breakpoints = rand_breakpoints(rng, recmap, x[i])
        x_[k] = recombination(x[i], x[N+i], breakpoints)  
        # second haplotype of indiv. k is the result of recombination of two
        # haplotypes of individual j (stored at indices j and N+j)
        breakpoints = rand_breakpoints(rng, recmap, x[j])
        x_[N+k] = recombination(x[j], x[N+j], breakpoints) 
        mutation!(rng, x_[k]  , arch)
        mutation!(rng, x_[N+k], arch)
    end
    DiploidWFPopulation(x_, x, N, arch, recmap) 
end

function remove_fixed!(pop::DiploidWFPopulation{<:NeutralArchitecture})
    idx = findall(x->x == 2pop.N, sum(pop.x))
    for g in pop.x
        g[idx] .= false
        dropzeros!(g)
    end
end

