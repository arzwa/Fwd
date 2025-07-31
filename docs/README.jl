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

# ## A purely neutral simulation
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


# ## Barrier locus
using Plots, WrightDistribution 

NA = 100
NB = 500
L = 1
C = 0.2
s = 0.05
m = 0.005
h = 0.5
u = s*h/200
xs = collect(C/2L:C/L:C)
AA = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
AB = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u  ) for _=1:L], xs)
R  = LinearMap(C)
xA = [ones(Bool, L) for _=1:2NA]
xB = [zeros(Bool, L) for _=1:2NB]
nA = collect(1:2NA)
nB = collect(1:2NB) .+ 2NA
popA = DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=deepcopy(xA), nodes=nA)
popB = DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=deepcopy(xB), nodes=nB)

qs, ts, pts = let
    rng  = Random.seed!(1)
    mpop = TwoPopOneWay(m, popA, popB)
    ts   = init_ts(mpop, C) 
    qs   = Vector{Float64}[]
    ngen = 20NB
    for i=1:ngen
        mpop = Fwd.generation!(rng, mpop, ts)
        push!(qs, mean(mpop.popB.x))
        if i % NB == 0
            mpop, ts = Fwd.simplify!(mpop, ts)
        end
    end
    rts = reverse_relabel(ts)
    pts = Fwd.to_tskit(rts)
    qs, ts, pts
end

stephist(hcat(qs...)', norm=true)
d = Wright(-2NB*s, 2NB*u, 2NB*(m + u), h)
plot!(0:0.001:1, p->pdf(d,1-p))

