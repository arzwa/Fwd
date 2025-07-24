using Fwd, Plots, Test, BenchmarkTools, WrightDistribution, Distributions
using Serialization, Barriers
default(grid=false, fontfamily="Computer modern", framestyle=:box,
    legend=false, fg_legend=:transparent)

# Parameters --------------------------------------------------------
L   = 100         # selected sites
G   = 10_000 + L  # sites
N   = 500         # population size
C   = 10.0        # Morgan
Ls  = 1.0 
s   = Ls/L
u   = s/200
m   = s
#dfe = Exponential(s)
dfe = Dirac(s)

# Seed the RNG
rng = Random.seed!(99)

# Random architecture -----------------------------------------------
xs = rand(G) .* C |> sort
ss = zeros(G)
selected = sort(sample(1:G, L, replace=false))
ss[selected] .= -rand(dfe, L)
qm = fill(1/2, G)
qm[selected] .= 1.0
neutral = map(i-> i ∉ selected, 1:G)
arch = Architecture([HaploidBiLocus(ss[i], u) for i=1:G], xs)

# Check scale of linkage between selected sites ---------------------
Rs = Fwd.rec_matrix(xs[selected]) 
heatmap(log2.(Rs ./ s), colorbar=true, clim=(-5,5), 
    size=(400,350), colormap=:RdBu)
# appears ok.

r̄ = map(i->mean(vcat(Rs[1:i-1,i], Rs[i+1:end,i]))/s, 1:L)

# Define the population model ---------------------------------------
haps    = [[j ∈ selected ? false : rand() < 0.5 for j=1:G] for _=1:N]
recmap  = LinearMap(C)
island  = HaploidWFPopulation(haps, arch, recmap)
mainland= HaploidHWLEPopulation(qm) 
metapop = MainlandIsland(island, mainland, m)

# Conduct a simulation ----------------------------------------------
Xs = simulate!(metapop, 510N, 
    callback=P->sum(P.island.x)/P.island.N, every=10)

Xs = Xs[10N+2:end]

X  = hcat(Xs...)

map(sample((1:G)[neutral], 16, replace=false)) do k
    stephist(X[:,k], bins=0:0.02:1, title="locus $k")
end |> x->plot(x..., size=(700,600))

# Store the 'data' --------------------------------------------------
data = (qm = qm, sim = Xs, ss = ss, xs = xs, selected = selected, neutral = neutral)
serialize("data/haploid-example.jls", data)

(qm, Xs, ss, xs, selected, neutral) = deserialize("data/haploid-example.jls")

# Prediction --------------------------------------------------------
ssel = ss[selected]
bloci = [Barriers.DiploidLocus(-2ssel[i], 0.5, u) for i=1:L]
barch = Barriers.Architecture(bloci, xs[selected], 
    Barriers.recrates(xs[selected]))
p̄ = zeros(L)
bmodel = Barriers.MainlandIslandModel(barch, m, p̄, N)
eq = Barriers.Equilibrium(bmodel)
eqq = 1 .- eq.Ep

# allele frequencies
qs = mean(Xs)[selected]
P1 = scatter(eqq, qs, size=(230,220), ms=1.5, color=:black, 
    xlabel="\$\\mathbb{E}[q]\$", ylabel="\$\\hat{q}\$")
plot!(x->x, color=2, xlim=(0,1), ylim=(0,1))

# neutral Fst
# recall: 
#   Fst = Var(p)/pq = (E[p^2] - p^2)/pq = 1 - E[pq]/pq
fn2 = var(Xs)[neutral] ./ 0.25
wins = Fwd.winstat(mean, 0.01, xs[neutral], fn2)
P2 = plot(map(x->(x, 1/(1+2N*(u .+ Barriers.me(eq, x)))), 0:0.01:C), 
    size=(600,200), color=:black, ylim=(0,1), label="prediction")
plot!(map(mean, first.(wins)), last.(wins), lw=1, label="simulation",
    ylabel="\$F_\\mathrm{ST}\$", margin=4Plots.mm, xlabel="map position",
    legend=:bottomright)
plot(P1, P2, layout=grid(1,2,widths=[0.25,0.75]), size=(700,200))

# Take a sample -----------------------------------------------------
# Say we pick a sample of two mainland and two island haplotypes.
rng = Random.seed!(15)
sample1 = [Fwd.rand_migrant(rng, metapop.mainland) for _=1:2]
sample2 = sample(rng, metapop.island.x, 2, replace=false)
smple = hcat(sample1..., sample2...)[neutral,:]

serialize("data/haploid-example-sample.jls", (xs=xs[neutral], ys=smple))


