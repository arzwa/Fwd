abstract type Locus end

mutationrate(l::Locus) = l.u

struct BiAllelic{T} <: Locus
    u :: T
end

struct IntAllelic{T} <: Locus
    u :: T
end

# XXX should we add the recombination map as a field?
struct Architecture{L,V<:AbstractVector,R<:RecombinationMap}
    loci   :: Vector{L}  # loci
    xs     :: V  # map locations
    recmap :: R
end
Base.length(arch::Architecture) = length(arch.loci)
Base.getindex(arch::Architecture, i) = arch.loci[i]

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

function mutation!(rng, locus::BiAllelic, X, N, i)
    u = mutationrate(locus)
    n = rand(rng, Poisson(u*N))
    idx = sample(rng, 1:N, n)  # haplotypes
    for j in idx
        X[j][i] = (X[j][i] == 1 ? 0 : 1)  # works for both Bool and Int
    end
end

function mutation!(rng, locus::IntAllelic, X, N, i)
    u = mutationrate(locus)
    ns = rand(rng, Poisson(u), N)
    for j=1:N
        X[j][i] += ns[j]
    end
end

