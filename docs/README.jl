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
# faster (?) compared to the `tskit` implementation), but of course it
# would be less prone to bugs...
#
# The library currently supports the simulation of Wright-Fisher
# populations
# - with linear genetic maps 
# - with multiple chromosomes 
# - embedded in a metapopulation
# Currently, there is only support for soft selection in the latter. The
# metapopulation implements migration as copying (at the beginning of each
# generation, before selection, an expected proportion mᵢⱼ of population j
# is replaced by clones send out by population i).

using Fwd, Random, StatsBase, Plots

# ## A purely neutral simulation

let N=3, C=0.1,   # popsize, map length, recombination map
    arch = Architecture(DiploidLocus{Float64}[], Float64[], LinearMap(C))
    pop = WFPopulation(ploidy=Diploid(), N=N, arch=arch)
    rng = Random.seed!(19)
    pop, ts = Fwd.simulate!(rng, pop, init_ts(pop), 3)
    pts = to_tskit(Fwd.reverse_relabel(ts))
    sts = simplify(ts, pop.nodes, keep_roots=true)
    print(draw_text(ts))
    print(draw_text(sts))
    print(pts.simplify(0:2N-1, keep_input_roots=true).draw_text())
end


# ## A single barrier locus
 
# We simulate a pair of haploid populations.
NA = 100
NB = 500
s = 0.05
m = 0.005
u = s/200
C = 0.1
AA = Architecture([BiAllelic(0.0)], Float64[C/2], LinearMap(C))
AB = Architecture([BiAllelic(  u)], Float64[C/2], LinearMap(C))
xA = [[true]  for _=1:NA]
xB = [[false] for _=1:NB]
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
Φ  = GPMap([HaploidLocus(-s, 1)])
popA = WFPopulation(ploidy=Haploid(), gpm=Φ, N=NA, arch=AA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), gpm=Φ, N=NB, arch=AB, x=xB, nodes=nB)
ngen = 20NB

mpop, ts, qs = let
    rng  = Random.seed!(1)
    mpop = TwoPopOneWay(m, popA, popB)
    mpop, ts, qs = simulate!(rng, mpop, init_ts(mpop), ngen, pop->mean(pop.popB.x)[1])
end

# Compare deleterious allele frequencies against theoretical prediction
# (diffusion theory)
using WrightDistribution
q = vec(hcat(qs...)')
d = Wright(-2NB*s, NB*u, NB*(m + u), 0.5)
stephist(q, norm=true, color=:gray, fill=true, fillalpha=0.2, label="simulation")
plot!(0:0.001:1, x->pdf(d,1-x), label="diffusion theory", xlabel="\$q\$", ylabel="density")
savefig("docs/pl1.png") #src
# ![](docs/pl1.png)

# get cross-population coalescence times
xs, _, _, tab = diffdiv(ts) 
plot(xs, tab, line=:steppost, color=:gray, fill=true, fillalpha=0.2, )
vline!(AB.xs, lw=2)
gff(x1, x2, s) = 1/(1+s/Fwd.recrate(abs(x1-x2))) 
plot!(x->1/(m*gff(x, C/2, s)), ylim=(0,ngen), color=:black, lw=2, 
    ls=:dash, xlabel="\$x\$", ylabel="\$T_{AB}\$")
savefig("docs/pl2.png") #src
# ![](docs/pl2.png)


# ## Multilocus cline

# $D$ demes along a 1D habitat in which $L$ alleles are divergently
# selected across the two halves.
D = 20
N = 50
m = 0.1
L = 10
C = 0.5
xs = C/2L:C/L:(C-C/2L)
s = 0.02
arch = Architecture([BiAllelic(s/200) for _=1:L], xs, LinearMap(C))
Ps = map(1:D) do d
    sgn = d > D÷2 ? 1 : -1
    x0 = d > D÷2 ? [ones(Bool, L) for _=1:N] : [zeros(Bool, L) for _=1:N]
    Φd = GPMap([HaploidLocus(sgn*s, i) for i=1:L])
    nd = collect((d*N+1):((d+1)*N))
    Pd = WFPopulation(ploidy=Haploid(), gpm=Φd, N=N, arch=arch, x=x0, nodes=nd)
end 
M = zeros(D,D)
for i=1:D
    i > 1 && (M[i,i-1] = m/2)
    i < D && (M[i,i+1] = m/2)
end
pop = MetaPop(Ps, M)
rng = Random.seed!(89)
pop, ts = simulate!(rng, pop, init_ts(pop), 20000);

# Graph the allele frequencies of teh selected loci at the end of the
# simulation and show the mean coalescence times between the leftmost and
# rightmost deme.
qs = permutedims(hcat([mean(pop[i].x) for i=1:D]...))
P1 = plot(qs, ylabel="\$q\$", xlabel="deme", marker=true, ms=2)
xx, _, _, tab = diffdiv(ts, 1, D)
P2 = plot(xx, tab, line=:steppost, xlabel="map position", ylabel="\$T_{1,$D}\$", color=:gray)
vline!(xs)
plot(P1, P2, size=(600,220), margin=3Plots.mm)

savefig("docs/pl3.png") #src
# ![](docs/pl3.png)
#
# The barrier effect is nicely demonstrated.
