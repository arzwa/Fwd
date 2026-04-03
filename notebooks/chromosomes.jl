# Multi-chromosome map
using Distributed#; addprocs(10)
@everywhere using Random, Fwd, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Parameters
import TreeSequences as TS
import MCMCChains: mcse

# Say we have 5 × 100 cM chromosomes
rng = Random.seed!(562)
C = 0.1
K = 3
chrs = [LinearMap(C) for _=1:K]
M = Chromosomes(chrs)

# Now have L=50 or so, and randomly scatter the loci across chromosomes
L   = 100
α   = 1.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
xs  = ys .* (C*K)

# Selection
Ls = 0.2
s̄  = Ls/L
ss = rand(Exponential(s̄), L)
u  = s̄/500
Ns = 5.

AA = Architecture([BiAllelic(0.) for _=1:L], xs, M)
AB = Architecture([BiAllelic(u ) for _=1:L], xs, M)
R  = Fwd.rec_matrix(AA)
r̄  = Fwd.rbar(R)

ms  = 1.0
m   = ms*s̄
NA  = 1 
NB  = ceil(Int, Ns/s̄)
MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
nA  = collect(1:NA)
nB  = collect(1:NB) .+ NA
xA  = [ ones(Int, L) for _=1:NA]
xB  = [zeros(Int, L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 1000000

using BarriersBP
BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, rm=M, Ne=NB, u=u))
AB = deepcopy(BP); AB.Ep .= 1.
ZSF = Equilibrium(ZSF24(m=m, s=2ss, h=fill(0.5,L), 
    xs=xs, rm=M, Ne=NB, u=u, p̄=zeros(L)))
nn = 1000
#plot!(range(0,C*K,nn), x->1/BarriersBP.me(AB, x), yscale=:log10)
P1 = plot(range(0,C*K,nn), x->1/BarriersBP.me(ZSF, x), yscale=:log10, color=:black, label="ZSF24")
plot!(range(0,C*K,nn), x->1/BarriersBP.me(BP, x), yscale=:log10, color=:firebrick, alpha=0.7, label="BP")
plot!(size=(900,200), legend=:bottomright)
vline!(cumsum(length.(M.maps)), color=:lightgray, alpha=0.8, lw=2, xlim=(0,Inf), label="")
P2 = scatter(BP.Ep, ZSF.Ep, xlabel="\$\\mathbb{E}[p]\$ (BP)", 
    ylabel="\$\\mathbb{E}[p]\$ (ZSF24)", color=:black, ms=2, size=(400,370))
plot(P1, P2, layout=grid(1,2,widths=[0.8,0.2]), margin=4Plots.mm)


ts = init_ts(mpop)
Fwd.generation!(rng, mpop, ts)

pop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), ngen, 
    pop->mean(pop.popB.x), simplify=100, every=100)

plot(mean(qs))

_ts = TS._add_grand_ancestor(ts)
xx, ta, tb, tab = TS.diffdiv(_ts)

plot(xx, tab, size=(900,200), yscale=:log10, 
    color=:gray, alpha=0.4, framestyle=:default, xlim=(-0.001,C*K))
title!(@sprintf("\$L=%2d, L\\bar{s}=%.2f, m/\\bar{s}=%.2f, \\bar{r}=%.2f, \\bar{r}/\\bar{s}=%.2f\$",
    L, Ls, ms, r̄, r̄/s̄))
vline!(cumsum(length.(M.maps)), color=:black)
plot!(range(0,C*K,nn), x->1/BarriersBP.me(BP, x), yscale=:log10, color=:firebrick, alpha=0.7, label="BP")
plot!(range(0,C*K,nn), x->1/BarriersBP.me(ZSF, x), yscale=:log10, color=:orange, alpha=0.7, label="BP")
plot!(xlabel="map position", ylabel="\$T\$", margin=5Plots.mm)
sticks!(twinx(), xs, ss, color=:firebrick, framestyle=:default, 
    xlim=(-0.001,C*K), lw=2, alpha=0.3, ylabel="\$s\$")


scatter(BP.Ep, 1 .- mean(qs))
plot!(x->x, color=:gray, alpha=0.4)


plot(xx, tab, size=(800,200), yscale=:log10, 
    color=:lightgray, framestyle=:default,
    xlim=(-0.01,Inf), margin=6Plots.mm)
title!(@sprintf("\$L=%2d, L\\bar{s}=%.2f, m/\\bar{s}=%.2f, \\bar{r}=%.2f, \\bar{r}/\\bar{s}=%.2f\$",
    L, Ls, ms, r̄, r̄/s̄))
plot!(range(extrema(xx)..., 1000), x->1/BarriersBP.me(BP, x), color=:teal)
plot!(range(extrema(xx)..., 1000), x->1/BarriersBP.me(AB, x), color=:orange)
vline!(cumsum(length.(M.maps)), color=:black)
xlabel!("map position (M)")
ylabel!("\$T\$")
sticks!(twinx(), xs, ss, color=:firebrick, framestyle=:default, lw=3, alpha=0.3,
    xlim=(-0.01,Inf), ylabel="\$s\\mathbb{E}[p]\$")

