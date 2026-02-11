using Fwd, Distributions, Random, Plots

# Architecture
rng = Random.seed!(255)
Ls  = 0.2
L   = 20
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
u   = s̄/200
mBC = s̄/10
mCB = s̄/5
Ns  = 20. 
NB  = ceil(Int64, Ns/s̄)
NA  = NB
NC  = 2NB
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
C   = 0.5
G   = C*100*10^6   # (1cM/Mb)
xs  = ceil.(Int64, ys .* G)
R   = LinearPhysMap(C=C, G=G)

# Ancestral population: this one starts from `00...0`.
# maybe let the ancestral evolve neutrally at all L loci
arch = Architecture([BiAllelic(u) for _=1:L], xs, R)
gpmA = GPMap([HaploidLocus(0.0, i) for i=1:L])
nA = collect(1:NA)
xA = [rand(Bool, L) for _=1:NA]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=arch, gpm=gpmA, x=xA, nodes=nA)

ngen = 10*max(NB,NC)
popA, ts, qs = Fwd.simulate!(popA, init_ts(popA), ngen, 
    pop->mean(pop.x), simplify=100, every=1)
plot(permutedims(hcat(qs...))[:,1:10])

pops = Fwd.splitpop(popA, ts, (NB, NC))

gpmB = GPMap([HaploidLocus( ss[i], i) for i=1:L])
gpmC = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
pops[1] = reconstruct(pops[1], gpm=gpmB)
pops[2] = reconstruct(pops[2], gpm=gpmC)

mpop = Fwd.MetaPop(pops, [0. mBC; mCB 0.])

ngen = 20*max(NB,NC)
mpop, ts, qs = Fwd.simulate!(mpop, ts, ngen, 
    pop->(mean(pop.P[1].x), mean(pop.P[2].x)), 
    simplify=100, every=ngen÷500)

PB = hcat(first.(qs)...) |> permutedims
PC = hcat(last.(qs)...) |> permutedims

map(1:20) do k
    plot(PB[:,k])
    plot!(PC[:,k])
end |> x->plot(x..., layout=(5,4), size=(700,600))

pb = mean(PB, dims=1) |> vec
pc = mean(PC, dims=1) |> vec
scatter(pb, pc, ylim=(0,1), xlim=(0,1))

k = 8
stephist( PB[:,k], bins=0:0.05:1.05)
stephist!(PC[:,k], bins=0:0.05:1.05)

ts = TS._add_grand_ancestor(ts)
xx, tb, tc, tbc = TS.diffdiv(ts, 1, 2)

plot(xx, tbc, size=(1000,200), )
sticks!(twinx(), xs, pb .- pc, lw=5, color=2)

plot!(xx, tb)
plot!(xx, tc)

pts = TS.to_tskit(ts)

pts.samples(population=0)

fst = 1 .- ((tb .+ tc) ./ 2) ./ tbc
plot(xx, fst, size=(1000,200), )
sticks!(twinx(), xs, pb .- pc, lw=5, color=2)

BP = Equilibrium(BPModel(m=mBC, xs=ys*C, Ne=NB, u=u, s=ss))

plot(xx, tbc, size=(1000,200), color=:gray, alpha=0.6)
mes = map(x->(G*x/C, 1/Barriers.me(BP,x)+NB), range(0,C,500))
plot!(mes, yscale=:log10, lw=2, color=:black)

# should do proper SC prediction but with mₑ 
