# XXX Currently assuming models without epistatic interactions, or interactions
# between genetics and life history/sex/...
# What is the elegant, functional way to implement interactions of that sort?
abstract type Locus end
abstract type MutationSampler end

mutationrate(l::Locus) = l.u

struct HaploidBiLocus{T} <: Locus
    s :: T  # selection coefficient
    u :: T  # mutation rate
end

fitnesseffect(l::HaploidBiLocus, x) = l.s*x

# general ploidy case (multiplicative effects)
fitnesseffect(l::HaploidBiLocus, args...) = sum(args) * l.s  

function mutation!(x::Vector{Bool}, i::Int) 
    x[i] = !x[i]
end

struct DiploidBiLocus{T} <: Locus
    s01 :: T
    s11 :: T
    u   :: T  
    # XXX this should be the diploid mutation rate with the current
    # logic of the mutation implementation.
end

fitnesseffect(l::DiploidBiLocus, x, y) = x != y ? l.s01 : (x == 1 ? l.s11 : 0.0)

# This is a model for an asexual stretch of genome, accumulating deleterious
# mutations with a multiplicative fitness effect. 
struct PoissonLocus{T}
    s :: T
    u :: T
end

fitnesseffect(l::PoissonLocus, x) = x*l.s
fitnesseffect(l::PoissonLocus, x, y) = (x+y)*l.s  # diploid case

struct Architecture{T,L,M<:MutationSampler}
    loci :: Vector{L}  # loci
    xs   :: Vector{T}   # map locations
    mut  :: M
end

Architecture(loci, xs) = Architecture(loci, xs, PoissonMutationSampler(loci))

logfitness(a::Architecture, x) = logfitness(a.loci, x)
function logfitness(a::Vector{L}, x) where L
    mapreduce(i->fitnesseffect(a[i], x[i]), +, 1:length(x)) 
end

fitness = exp ∘ logfitness

struct PoissonMutationSampler{P<:Poisson,W<:Weights} <: MutationSampler
    d  :: P
    ws :: W
end

function PoissonMutationSampler(loci::Vector{L}) where L<:Locus
    us = mutationrate.(loci) 
    PoissonMutationSampler(Poisson(sum(us)), Weights(us))
end

Base.length(M::PoissonMutationSampler) = length(M.ws)

function rand_mutations(rng::AbstractRNG, M::PoissonMutationSampler)
    nmut = rand(rng, M.d)
    mutations = sample(rng, 1:length(M), Weights(M.ws), nmut)
end

function mutation!(x::Vector, mutations::Vector{Int})
    for i in mutations
        mutation!(x, i)
    end
    return x
end




