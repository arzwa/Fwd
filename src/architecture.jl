abstract type Locus end
abstract type MutationSampler end

mutationrate(l::Locus) = l.u

struct HaploidBiLocus{T} <: Locus
    s :: T  # selection coefficient
    u :: T  # mutation rate
end

fitnesseffect(l::HaploidBiLocus, x) = l.s*x

struct DiploidBiLocus{T} <: Locus
    s01 :: T
    s11 :: T
    u   :: T  
end

fitnesseffect(l::DiploidBiLocus, x, y) = x != y ? l.s01 : (x == 1 ? l.s11 : 0.0)

# This is a model for an asexual stretch of genome, accumulating deleterious
# mutations with a multiplicative fitness effect. 
@with_kw struct PoissonLocus{T} <: Locus
    @assert s > 0.0
    s :: T 
    u :: T
end

fitnesseffect(l::PoissonLocus, x) = -x*l.s
fitnesseffect(l::PoissonLocus, x, y) = -(x+y)*l.s  # diploid case

struct Architecture{L,V<:AbstractVector}#M<:MutationSampler}
    loci :: Vector{L}  # loci
    xs   :: V  # map locations
end
Base.length(arch::Architecture) = length(arch.loci)
Base.getindex(arch::Architecture, i) = arch.loci[i]

Architecture() = Architecture(HaploidBiLocus{Float64}[], Float64[])

logfitness(a::Architecture, args...) = logfitness(a.loci, args...)

function logfitness(a::Vector{L}, x) where L
    length(a) == 0 && return 0.0
    mapreduce(i->fitnesseffect(a[i], x[i]), +, 1:length(x)) 
end

function logfitness(a::Vector{L}, x::Tuple) where L
    length(a) == 0 && return 0.0
    mapreduce(i->fitnesseffect(a[i], x[1][i], x[2][i]), +, 1:length(x)) 
end

fitness = exp ∘ logfitness

# NOTE: Mutation is implemented at the population level, i.e. the outer
# loop is over loci, for each locus we apply mutation to the entire
# population.  Earlier, we had mutation at the individual level, which
# seems nicer from an individual-based model point of view, but in fact
# isn't when the `Architecture` is defined at the population level!
# Arguably, the right way to go is to enable arbitrary subpopulations with
# different architectures (this is straightforward to implement I think
# when we would happen to need it). In the latter case, in a `generation`
# function one would then simply iterate over `(X,A)` pairs and apply
# `mutation!` to those.
function mutation!(rng, X, A::Architecture)
    N = length(X)
    for i=1:length(A)
        mutation!(rng, A[i], X, N, i)
    end
    return X
end

function mutation!(rng, locus::Union{HaploidBiLocus,DiploidBiLocus}, X, N, i)
    u = mutationrate(locus)
    n = rand(rng, Poisson(u*N))
    idx = sample(rng, 1:N, n)  # haplotypes
    for j in idx
        X[j][i] = (X[j][i] == 1 ? 0 : 1)  # works for both Bool and Int
    end
end

function mutation!(rng, locus::PoissonLocus, X, N, i)
    u = mutationrate(locus)
    ns = rand(rng, Poisson(u), N)
    for j=1:N
        X[j][i] += ns[j]
    end
end
