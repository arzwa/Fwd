# Tree sequence data structure
struct TreeSequence{T}
    nodes :: Vector{Int64}
    edges :: Vector{Edge}
    sites :: Vector{Site{T}}
    mutations :: Vector{Mutation{T}}
end

struct Edge
    left   :: Int
    right  :: Int
    parent :: Int
    child  :: Int
end

struct Site{T}
    position :: Int
    ancestral :: T
end

struct Mutation{T}
    site :: Int
    node :: Int
    derived :: T
end

# Forward simulation
abstract type Locus end

struct HaploidLocus{T} <: Locus
    s :: T  # selection coefficient
    u :: T  # mutation rate
end

# Genomic architecture
struct Architecture{T<:Real,V<:Locus,R<:RecombinationMap}
    U      :: T            # total mutation rate
    us     :: Vector{T}    # mutation rates of modeled sites
    loci   :: Vector{V}    # vector of loci
    recmap :: R            # XXX recombination.jl
end


