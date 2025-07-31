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
export TreeSequence, reverse_relabel, simplify, to_tskit, from_tskit, draw_text

include("diploids.jl")
export DiploidWFPopulation, generation!, TwoPopOneWay, init_ts

include("utils.jl")

end # module Fwd


