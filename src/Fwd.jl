module Fwd

using SparseArrays, Random, Reexport, Distributions, Parameters
using LinearAlgebra, StatsBase
using ProgressMeter
import Random: AbstractRNG
@reexport using Random
@reexport using SparseArrays
@reexport using Random: default_rng

const Gb = 1_000_000_000
const Mb = 1_000_000
const kb = 1_000

export Gb, Mb, kb

include("architecture.jl")
include("haploids.jl")
include("mainland.jl")

end # module Fwd
