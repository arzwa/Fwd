# Two (explicit) population model, with unidirectional migration, **and
# selection**.
using Pkg; Pkg.activate("/home/arzwa/dev/Fwd")
using Fwd, Plots, Test, BenchmarkTools, WrightDistribution, Distributions
using Serialization, Barriers, Parameters
plotsdefault()

# Fst function for two-island/unidirectional migration --------------
fstfun = eval(
    :(function (m, N_A, N_B)
      #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =# @inbounds begin
              #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =#
              begin
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:389 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:390 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:391 =#
                  (/)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(1//2, N_A), m)), (*)((*)(1//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-1//1, N_A), (^)(m, 2))), (*)((*)(-1//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)(1//2, N_A), (^)(m, 3))), (*)((*)(1//2, N_B), (^)(m, 3))), (+)((+)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(3//2, N_A), m)), (*)((*)(7//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-3//1, N_A), (^)(m, 2))), (*)((*)(-5//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)((*)(4//1, N_A), N_B), (^)(m, 2))), (*)((*)(3//2, N_A), (^)(m, 3))), (*)((*)(3//2, N_B), (^)(m, 3))), (*)((*)((*)(-2//1, N_A), N_B), (^)(m, 3))))
              end
          end
  end)
)

# Parameters --------------------------------------------------------
L   = 100         # selected sites
G   = 10_000 + L  # sites
NA  = 500        # population size
NB  = 500
C   = 10.0        # Morgan
Ls  = 1.0 
s   = Ls/L
u   = s/500
m   = 0.8s
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
archA = Architecture([HaploidBiLocus(-ss[i], u) for i=1:G], xs)
archB = Architecture([HaploidBiLocus( ss[i], u) for i=1:G], xs)

# Check scale of linkage between selected sites ---------------------
Rs = Fwd.rec_matrix(xs[selected]) 
heatmap(log2.(Rs ./ s), colorbar=true, clim=(-5,5), 
    size=(400,350), colormap=:RdBu)
# appears ok.

r̄ = map(i->mean(vcat(Rs[1:i-1,i], Rs[i+1:end,i]))/s, 1:L)

# Define the population model ---------------------------------------
rng = Random.seed!(126)
hapsA   = [[j ∈ selected ? true  : rand() < 0.5 for j=1:G] for _=1:NA]
hapsB   = [[j ∈ selected ? false : rand() < 0.5 for j=1:G] for _=1:NB]
recmap  = LinearMap(C)
island  = HaploidWFPopulation(hapsB, archB, recmap)
mainland= HaploidWFPopulation(hapsA, archA, recmap) 
metapop = MainlandIsland(island, mainland, m)

# Conduct a simulation ----------------------------------------------
Xs = simulate!(metapop, 1100NA, 
    callback=P->(allele_freqs(P.mainland), allele_freqs(P.island)), 
    every=25)

Xi = hcat(last.(Xs)[end-10000+1:end]...)
Xm = hcat(first.(Xs)[end-10000+1:end]...)

map(sample((1:G)[neutral], 16, replace=false)) do k
    stephist(Xm[k,:], bins=0:0.02:1.02, title="locus $k")
    stephist!(Xi[k,:], bins=0:0.02:1.02, title="locus $k")
end |> x->plot(x..., size=(700,600))

# Store the 'data' --------------------------------------------------
fname = "haploid-example-2b-G$G"
data = (pop=metapop, qm=qm, sim=Xs, 
    NA=NA, NB=NB, m=m, u=u, ss=ss, xs=xs, 
    selected=selected, neutral=neutral)
serialize("data/$fname.jls", data)

# load the data -----------------------------------------------------
fname = "haploid-example-2b-G$G"
data = deserialize("data/$fname.jls")
@unpack ss, neutral, selected, xs, m, u, NA, NB, sim = data
Xi = hcat(last.(sim)[end-10000+1:end]...)
Xm = hcat(first.(sim)[end-10000+1:end]...)

# Prediction --------------------------------------------------------
ssel = ss[selected]
bloci = [Barriers.DiploidLocus(-2ssel[i], 0.5, u) for i=1:L]
barch = Barriers.Architecture(bloci, xs[selected], 
    Barriers.recrates(xs[selected]))
p̄ = zeros(L)
bmodel = Barriers.MainlandIslandModel(barch, m, p̄, NB)
eq = Barriers.Equilibrium(bmodel)
eqq = 1 .- eq.Ep

# allele frequencies
qs = mean(Xi, dims=2)[selected]
P1 = scatter(1 .- eqq, 1 .- qs, size=(230,220), ms=1.5, color=:black, 
    xlabel="\$\\mathbb{E}[\\Delta]\$", ylabel="\$\\hat{\\Delta}\$", legend=false)
plot!(x->x, color=2, xlim=(0,1), ylim=(0,1))

# Pair Fst
pii = 2Xi .* (1 .- Xi)
pim = 2Xm .* (1 .- Xm)
pib = (Xi .* (1 .- Xm) .+ (1 .- Xi) .* Xm) 
piw = (pii .+ pim) ./ 2

# should one do this:
#fim = (pib .- piw) ./ (pib .+ piw)
#fst = map(eachrow(fim)) do row
#    mean(filter(!isnan, row))
#end
# or rather this:
pib_ = mean(pib, dims=2) 
piw_ = mean(piw, dims=2) 
fst = (pib_ .- piw_) ./ (pib_ .+ piw_)
fst = fst[neutral,:]
wins = Fwd.winstat(mean, 0.05, xs[neutral], fst)
# I guess for the sake of comparing to theoretical predictions based on
# expected coalescence times, the latter is better.

# Is this what people do in practice? Does it in expectation yield the fst
# estimate above?
fstgens = map(1:1:size(Xi, 2)) do gen
    piiw = Fwd.winstat(mean, 0.05, xs[neutral], piw[neutral,gen])
    piib = Fwd.winstat(mean, 0.05, xs[neutral], pib[neutral,gen])
    a = last.(piiw); b = last.(piib)
    fstw = (b .- a) ./ (a .+ b)
end

P2 = plot(map(x->(x, fstfun(Barriers.me(eq, x), NA, NB)), 0:0.01:C), 
    size=(600,200), color=:black, ylim=(0,1), label="prediction", xlim=(0,C))
vline!(xs[selected], alpha=0.2, color=:gray, label="")
plot!(map(mean, first.(wins)), last.(wins), lw=1, label="simulation",
    ylabel="\$F_\\mathrm{ST}\$", margin=4Plots.mm, xlabel="map position",
    legend=:topright)
#plot!(map(mean, first.(wins)), mean(fstgens), ms=1, alpha=0.3,
#    label="simulation, windowed", ylabel="\$F_\\mathrm{ST}\$", margin=4Plots.mm, 
#    xlabel="map position", legend=:topright)
plot!(map(mean, first.(wins)), fstgens[end], ms=1, alpha=0.7,
    label="last generation, windowed", legend=:outertopright)
P0 = plot(P1, P2, layout=grid(1,2,widths=[0.10,0.90]), size=(1100,220), 
    margin=6Plots.mm)

# Coarse predictions ---------------------------------------------
using StatsBase
nwin = 200
Δ = step(range(0, C, nwin))
X = fit(Histogram, xs[selected], 0:Δ:C+Δ).weights
CM = Barriers.CoarseModel(X=X, Δ=fill(Δ, nwin), s=s, m=m, u=u, λ=1/NA)
mec = m .* Barriers.gff(CM)
CM = Barriers.CoarseModel(X=X, Δ=fill(Δ, nwin), s=s, m=m, u=u, λ=0.0)
mec2 = m .* Barriers.gff(CM)
P3 = plot(0:0.01:C, map(x->Barriers.me(eq, x), 0:0.01:C), xlim=(0,C),
    xlabel="map position (M)", margin=6Plots.mm, ylabel="\$m_e\$",
    size=(1000,200), color=:black, label="prediction", legend=:outertopright)
vline!(xs[selected], alpha=0.2, color=:gray, label="")
plot!(0:Δ:C+Δ, [mec ; mec[end]], linetype=:steppost, lw=1, label="coarse, partial div.")
plot!(0:Δ:C+Δ, [mec2 ; mec2[end]], linetype=:steppost, lw=1, label="coarse, complete div.")

# it seems we did not get the divergence right, but almost so -- this does
# however make quite a difference.

P4 = plot(map(x->(x, fstfun(Barriers.me(eq, x), NA, NB)), 0:0.01:C), 
    size=(600,200), color=:black, ylim=(0.2,0.8), label="prediction", xlim=(0,C),
    ylabel="\$F_\\mathrm{ST}\$", margin=4Plots.mm, xlabel="map position",
    legend=:outertopright)
plot!(0:Δ:C+Δ, fstfun.([mec ; mec[end]], NA, NB), 
    linetype=:steppost, lw=1, label="coarse, partial div.")
plot!(0:Δ:C+Δ, fstfun.([mec2 ; mec2[end]], NA, NB), 
    linetype=:steppost, lw=1, label="coarse, complete div.")

plot(P0, P3, P4, layout=(3,1), size=(1100, 600))

# Take a sample -----------------------------------------------------
# Say we pick a sample of `n` mainland and `n` island haplotypes.
rng = Random.seed!(15)
n = 10
sample1 = sample(rng, metapop.mainland.x, n, replace=false)
sample2 = sample(rng, metapop.island.x, n, replace=false)
smple = hcat(sample1..., sample2...)[neutral,:]

serialize("data/$fname-sample.jls", (xs=xs[neutral], ys=smple))

# just store the entire population for later use
rng = Random.seed!(156)
mm = hcat(metapop.mainland.x...)[:,shuffle(1:NA)]
ii = hcat(metapop.island.x...)[:,shuffle(1:NB)]
pops = [mm ii][neutral,:]

popdata = (xs=xs[neutral], X=pops, NA=NA, NB=NB)

serialize("data/$fname-pop.jls", popdata)
