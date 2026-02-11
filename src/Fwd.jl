__precompile__()
module Fwd

using Random
import Random: AbstractRNG
using Distributions
using Parameters
#using LinearAlgebra
using StatsBase
using Printf
using ProgressMeter
using TreeSequences
import TreeSequences as TS

const Gb = 1_000_000_000
const Mb = 1_000_000
const kb = 1_000
export Gb, Mb, kb

include("rec.jl")
export LinearMap, LinearPhysMap, Unlinked, Chromosomes
export maplength, rand_breakpoints

include("arch.jl")
include("gpm.jl")
export BiAllelic, IntAllelic, HaploidLocus, HaploidTwoLocus, DiploidLocus
export Architecture, GPMap

#include("ts.jl")
#export TreeSequence, reverse_relabel, simplify, to_tskit, from_tskit, draw_text

include("wfpop.jl")
export WFPopulation, Haploid, Diploid, generation!, init_ts

include("twopop.jl")
export TwoPopOneWay

include("metapop.jl")
export MetaPop

include("utils.jl")
export diffdiv, simulate!

end # module Fwd


