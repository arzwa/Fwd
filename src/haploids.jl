abstract type Population end

struct HaploidWFPopulation{A<:Architecture,R<:RecombinationMap,T} <: Population
    x  :: Vector{T}
    x_ :: Vector{T}
    N  :: Int64
    arch :: A
    recmap :: R
end
    
function HaploidWFPopulation(x, arch, recmap)
    HaploidWFPopulation(x, similar(x), length(x), arch, recmap)
end

function count_alleles(xs::Vector{SparseVector{Bool,Int}})
    d = Dict{Int64,Int64}()
    for x in xs
        for i in x.nzind
            if !haskey(d, i)
                d[i] = 1
            else
                d[i] += 1
            end
        end
    end
    ks = collect(keys(d))
    vs = collect(values(d))
    o = sortperm(ks)
    SparseVector(length(xs[1]), ks[o], vs[o])
end

#allelefreqs(pop::HaploidWFPopulation) = sum(pop.x) ./ pop.N
allelefreqs(pop::HaploidWFPopulation) = count_alleles(pop.x) ./ pop.N

function Base.show(io::IO, pop::HaploidWFPopulation)
    @unpack x, N, arch, recmap = pop
    write(io, "HaploidWFPopulation\n  $arch\n  $recmap\n  " *
        "L = $(length(x[1])/Mb)Mb\n  N = $N")
end

function simulate(rng::AbstractRNG, pop, n, args...; 
        cb=allelefreqs, show_progress=true, remove_fixed=Inf)
    x0 = cb(pop)
    xs = Array{typeof(x0)}(undef, n); xs[1] = x0
    pr = Progress(n; enabled=show_progress)
    for i=2:n
        pop = generation!(rng, pop, args...)
        xs[i] = cb(pop)
        i % remove_fixed == 0 && remove_fixed!(pop)
        next!(pr)
    end
    return xs
end

function generation!(rng::AbstractRNG, 
        pop::HaploidWFPopulation)
    @unpack N, x, x_, arch, recmap = pop
    w   = map(xi->fitness(arch, xi), x)
    idx = sample(rng, 1:N, Weights(w), 2N)
    @inbounds for k=1:N
        i = idx[k]    # parent 1 index
        j = idx[N+k]  # parent 2 index
        # first haplotype of indiv. k is the result of recombination of two
        # haplotypes of individual i (stored at indices i and N+i)
        breakpoints = rand_breakpoints(rng, recmap, x[i])
        x_[k] = recombination(x[i], x[j], breakpoints)  
        mutation!(rng, x_[k], arch)
    end
    HaploidWFPopulation(x_, x, N, arch, recmap) 
end

function remove_fixed!(pop::HaploidWFPopulation)
    @unpack arch = pop
    for (region, l) in zip(arch.xs, arch.loci)
        !isneutral(l) && continue
        remove_fixed_region!(pop, region)
    end
end

function remove_fixed_region!(pop, region::Region{V}) where V<:AbstractVector{Int}
    cts = count_alleles(pop.x)
    idx = findall(x->x == pop.N, cts[region.xs])
    for genome in pop.x
        genome[idx] .= false
        dropzeros!(genome)
    end
end

function remove_fixed_region!(pop, region::Region) 
    cts = count_alleles(pop.x)
    idx = vcat([findall(x->x == pop.N, cts[subregion]) 
        for subregion in region.xs]...)
    for genome in pop.x
        genome[idx] .= false
        dropzeros!(genome)
    end
end

