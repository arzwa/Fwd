# # Fwd
# 
# Forward population genetics simulation in julia, with tree sequence
# recording. Motivated by not wanting to force SLiM into simulating the
# sorts of models it was not designed for.
#
# The package implements a basic tree sequence data structure, and
# functions to convert to and from `tskit` tree sequences. `tskit` and
# `msprime` can then be used for recapitation and downstream analyses.
#
# I guess one could work directly with the `tskit` data structures using
# the low level implementations in the latter. Speed-wise that should not
# yield much of a difference (benchmarks suggest tree sequence
# simplification as implemented here is about equally fast if not slightly
# faster compared to the `tskit` implementation), but of course it would be
# less prone to bugs...

# Example of a purely neutral simulation
using Fwd, Random, StatsBase

# Set up:
N = 3  # diploid individuals
C = 0.1  # map length
A = Architecture(DiploidBiLocus{Float64}[], Float64[])  # empty selective architecture
R = LinearMap(C)
x = [Bool[] for _=1:2N] # haplotypes for selected loci
ts = Fwd.init_ts(2N, C) 

let 
    rng = Random.seed!(19)
    pop = DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    for _=1:3
        idx = sample(rng, 1:N, 2N)
        pop = Fwd.generation!(rng, pop, idx, ts)
    end
    pts = to_tskit(Fwd.reverse_relabel(ts))
    sts = simplify(ts, pop.nodes, keep_roots=true)
    print(draw_text(ts))
    print(draw_text(sts))
    print(pts.simplify(0:2N-1, keep_input_roots=true).draw_text())
end
