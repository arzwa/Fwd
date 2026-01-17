
using Distributed;  addprocs(5)
@everywhere using Random, Fwd, ProgressMeter, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Barriers

rng = Random.seed!(15)
Ls  = 0.4
L   = 40
s   = Ls/L
u   = s/200
ms  = 2
m   = ms*s
rs  = 1.0
r   = rs*s
xs  = cumsum([Fwd.distance(r) for i=1:L])
C   = xs[end] + xs[1] 
NA  = 1
Ns  = 5. 
NB  = ceil(Int64, Ns/s)
loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
EM = Equilibrium(Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1))
AM = AeschbacherModel(m, fill(s,L), xs)
BP = Equilibrium(BPModel(m, fill(s,L), xs, NB, u))
plot(range(0, C, 500), x->Barriers.me(EM, x))
plot!(range(0, C, 500), x->Barriers.me(AM, x))
plot!(range(0, C, 500), x->Barriers.me(BP, x))

plot(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10)
plot!(range(0, C, 500), x->1/Barriers.me(AM, x))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x))

# take ngen ~ 10 × max cross-pop coalescence time
minme = quantile(map(x->Barriers.me(BP, x), range(0, C, 100)), 0.05)
ngen = ceil(Int, (10 / minme * 1e-4)) * 10^4  

ngen = 500_000
nrep = 5

res = pmap(1:nrep) do _
    R   = LinearMap(C)
    AA  = Architecture([BiAllelic(0.0)   for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u)     for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
    nA  = collect(1:NA)
    nB  = collect(1:NB) .+ NA
    xA  = [ ones(Int, L) for _=1:NA]
    xB  = [zeros(Int, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    pop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), ngen, 
        pop->mean(pop.popB.x), simplify=100, every=100)
end

PP = map(x->1 .- permutedims(hcat(x...)), last.(res))

ps = vec(mean(mean(PP), dims=1))
P0 = scatter(ps, EM.Ep, color=:black, ms=2, xlim=(0.5,1), ylim=(0.5,1))
plot!(x->x)
P1 = scatter(ps, BP.Ep, color=:black, ms=2, xlim=(0.5,1), ylim=(0.5,1))
plot!(x->x)
plot(P0, P1, size=(460,220), margin=3Plots.mm, 
    title="\$r/s=$rs,\ m=$m "
    ylabel="\$\\hat{p}\$", xlabel="\$\\bar{p}\$")


