using Fwd, Plots, Test, BenchmarkTools, WrightDistribution
using Barriers, StatsBase, LinearAlgebra
default(grid=false, fontfamily="Computer modern", framestyle=:box, legend=false)

# Himani's example from fig 1 (Sachdeva 2022)
N  = 100
C  = 500.
s  = 0.02
u  = s*0.005
m  = s*0.45
L  = 80
xs = collect(range(0, C, L))
R  = Fwd.rec_matrix(xs)
M  = LinearMap(C)
haps    = [[rand() < u/s for i=1:L] for _=1:N]
arch    = Architecture([HaploidBiLocus(-s, u) for i=1:L], xs)
island  = HaploidWFPopulation(haps, arch, M)
mainland= HaploidFixedPopulation(fill(true, L)) 
metapop = MainlandIsland(island, mainland, m)

Xs = simulate!(metapop, 100N)

Xs = simulate!(metapop, 5100N, callback=P->sum(P.island.x)/P.island.N)
X = hcat(Xs...)[:,100N+2:(N÷2):end]

bloci  = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
barch  = Barriers.Architecture(bloci, xs, Barriers.recrates(xs))
p̄      = zeros(L)
bmodel = Barriers.MainlandIslandModel(barch, m, p̄, N)
eq     = Barriers.Equilibrium(bmodel)
d      = Barriers.distributions(eq)

function loghist(X, bins=0:0.01:1.01)
    hist = normalize(fit(Histogram, vec(X), bins))
    es = hist.edges[1]
    ws = hist.weights
    ws[end-1] += ws[end]
    ys = log.(ws[1:end-1])
    xs = [es[i] + step(es)/2 for i=1:length(ys)-1]
    return xs, reverse(ys)
end

scatter(loghist(X, 0:0.01:1.01)..., ms=2, ylim=(-Inf,Inf), label="")
plot!(0:0.0001:1, p->logpdf(d[1], p)) 


# -------------------------------------------
ms = 0.15:0.1:0.75
res = map(ms) do m
    @info m*s
    L = 80
    xs = collect(range(0, C, L))
    R  = Fwd.rec_matrix(xs)
    rm = LinearMap(C)
    rng     = default_rng()
    haps    = [[rand() < u/s for i=1:L] for _=1:N]
    arch    = Architecture([HaploidBiLocus(-s, u) for i=1:L], xs)
    island  = HaploidWFPopulation(haps, arch, rm)
    mainland= HaploidFixedPopulation(fill(true, L)) 
    metapop = MainlandIsland(island, mainland, m*s)
    Xs = simulate!(metapop, 5100N, callback=P->sum(P.island.x)/P.island.N)
    hcat(Xs...)[:,100N+2:end]
end

sims = [(0.1, 0.9438), (0.15, 0.9265), (0.2, 0.9134), (0.25, 0.8917), (0.3, 0.872), (0.35, 0.8243), (0.4, 0.8076), (0.45, 0.7307), (0.5, 0.1272), (0.55, 0.0788), (0.6, 0.0633), (0.65, 0.0436), (0.7, 0.0343), (0.75, 0.028), (0.8, 0.024)]

bloci  = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
barch  = Barriers.Architecture(bloci, xs, Barriers.recrates(xs))
p̄      = zeros(L)
pred = map(minimum(ms):0.01:maximum(ms)) do m
    bmodel = Barriers.MainlandIslandModel(barch, m*s, p̄, N)
    eq     = Barriers.Equilibrium(bmodel)
    m, eq.Ep[1]
end

scatter(sims, ms=2, color=:black)
plot!(pred, color=:gray, size=(300,250), xlabel="\$m/s\$",
    ylabel="\$\\mathbb{E}[p]\$", title="\$L=$L, s=$s, u=s/200\$")

# --------------------------------------------------
