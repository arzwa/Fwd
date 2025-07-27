__precompile__()
module Fwd

using Random, Reexport, Distributions, Parameters
import Random: AbstractRNG
using LinearAlgebra, StatsBase, Printf
using ProgressMeter, DataStructures
using PyCall
const tskit = PyNULL()

function __init__()
    copy!(tskit, pyimport("tskit"))
end

const Gb = 1_000_000_000
const Mb = 1_000_000
const kb = 1_000


export Gb, Mb, kb

include("architecture.jl")
export HaploidBiLocus, DiploidBiLocus, Architecture, fitness, logfitness

include("recombination.jl")
export LinearMap, maplength, rand_breakpoints

include("ts.jl")

#include("haploids.jl")
#export HaploidWFPopulation, HaploidFixedPopulation, MainlandIsland
#export HaploidHWLEPopulation
#export generation!, allele_freqs
#
#include("simulation.jl")
#export simulate!

end # module Fwd


# Notes
#
# 1. We focus on models with a discrete set of loci. A population should be a
# matrix with haplotypes in rows, loci in columns, and an index keeping which
# individual has which haplotypes.
