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

using Fwd, Random, StatsBase, Plots

# ## A purely neutral simulation
#
# We have to provide an 'empty' genetic architecture to do a neutral forward
# similation with tree sequence recording.

# Set up:
N = 3  # diploid individuals
C = 0.1  # map length
A = Architecture(DiploidBiLocus{Float64}[], Float64[])  # empty selective architecture
R = LinearMap(C)
x = [Bool[] for _=1:2N] # haplotypes for selected loci

let 
    pop = WFPopulation(ploidy=Diploid(), N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    ts  = Fwd.init_ts(pop, C) 
    rng = Random.seed!(19)
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


# ## Barrier locus/loci
# 
# We'll simulate a haploid pair of populations.
NA = 100
NB = 500
L = 1
C = 0.2
s = 0.05
m = 0.005
u = s/200
xs = collect(C/2L:C/L:C)
AA = Architecture([Fwd.HaploidBiLocus(0., 0.) for _=1:L], xs)
AB = Architecture([Fwd.HaploidBiLocus(-s,  u) for _=1:L], xs)
R  = LinearMap(C)
xA = [ones( Bool, L) for _=1:NA]
xB = [zeros(Bool, L) for _=1:NB]
nA = collect(1:NA)
nB = collect(1:NB) .+ 2NA
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, recmap=R, x=deepcopy(xA), nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, recmap=R, x=deepcopy(xB), nodes=nB)
ngen = 20NB

qs, ts, pts = let
    rng  = Random.seed!(1)
    mpop = TwoPopOneWay(m, popA, popB)
    ts   = init_ts(mpop, C) 
    qs   = Vector{Float64}[]
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

# Compare deleterious allele frequencies against theoretical prediction
# (diffusion theory)
using WrightDistribution
q = vec(hcat(qs...)')
d = Wright(-2NB*s, NB*u, NB*(m + u), 0.5)
stephist(q, norm=true, color=:gray, fill=true, fillalpha=0.2, label="simulation")
plot!(0:0.001:1, x->pdf(d,1-x), label="diffusion theory", xlabel="\$q\$", ylabel="density")
savefig("data/pl1.svg") ##
# ![](data/pl1.svg)

# get tree heights
pts = pts.simplify()
heights = map(tree->[tree.time(r) for r in tree.roots], pts.trees())[1:end-1]
heights = map(h->length(h) > 1 ? ngen : h[1], heights)
bps = collect(pts.breakpoints())[2:end-1]

plot(bps, heights, linetype=:steppre, color=:black, xlabel="map position", 
    ylabel="tree height")
savefig("data/pl2.svg") ##
# ![](data/pl2.svg)
