

using Distributed #; addprocs(10)
@everywhere using Random, Fwd, ProgressMeter, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Barriers

rng = Random.seed!(672)
Ls  = 0.3
L   = 100
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)

C   = 0.5
α   = 1.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] .* C
xs  = [(zs[i] + zs[i+1])/2 for i=1:L]

u    = s̄/200
m    = s̄
NA   = 1
NB   = 500
    
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
EM = Barriers.Equilibrium(M)
AM = AeschbacherModel(m, [-ss[i] for i=1:L], xs)
plot(range(0, C, 500), x->Barriers.me(EM, x))
plot!(range(0, C, 500), x->Barriers.me(AM, x))

# take ngen ~ 10 × max cross-pop coalescence time
minme = minimum(map(x->Barriers.me(EM, x), range(0, C, 1000)))
ngen = ceil(Int, (10 / minme * 1e-4)) * 10^4  
nrep = 50

res = pmap(1:nrep) do _
    AA  = Architecture([BiAllelic(0.0)   for _=1:L], xs)
    AB  = Architecture([BiAllelic(u)     for _=1:L], xs)
    MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
    MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
    R   = LinearMap(C)
    nA  = collect(1:NA)
    nB  = collect(1:NB) .+ NA
    xA  = [ ones(Int, L) for _=1:NA]
    xB  = [zeros(Int, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, 
        gpm=MA, recmap=R, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, 
        gpm=MB, recmap=R, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    seed = rand(rng, 1:2^32)
    pop, ts, qs = Fwd.simulation!(seed, mpop, C, L, ngen)
    (seed, ts, qs)
end

xy = map(X->Fwd.diffdiv(X[2])[[1,4]], res)

Q = mapreduce(x->vec(mean(x, dims=1)), hcat,getindex.(res, 3))

xy, Q = deserialize("data/divsel-polygenic-2.jls")
x, Y = Fwd.summarize_wins(first.(xy), last.(xy))
y = vec(mean(Y, dims=1))

P1 = plot(x, 1 ./ (y .- NA), linetype=:steppre, color=:lightgray, alpha=0.5)
plot!(range(0, C, 500), x->Barriers.me(EM, x), color=:black)
plot!(range(0, C, 500), x->Barriers.me(AM, x))
plot!(size=(700,200),  
    ylabel="\$m_e\$",
    title="\$L\\bar{s} = $Ls, L=$L, N_e=$NB, m=$m\$")
P2 = sticks(xs, ss, color=:firebrick, ylabel="\$s\$",
    xlabel="map position (M)")
plot(P1, P2, bottom_margin=4Plots.mm, 
    layout=grid(2,1,heights=[0.7,0.3]), size=(700,300))

# These approaches are not working. Somehow, to get in the right me range,
# substituting the mean `s`, weighted by predicted divergence, for s in the
# original coarse approximation works...
nwin = 200
Δ = step(range(0, C, nwin+1))
Δs= fill(Δ, nwin)
X = fit(Histogram, xs, 0:Δ:C).weights

sd = mean(ss .* EM.Ep)

CM = Barriers.CoarseModel(X=X, Δ=Δs, s=sd, m=m, u=u, λ=0.0)
Barriers.predict_divergence(CM, 1/NB)
mec = m .* Barriers.gff(CM)
plot(range(0, C, 500), x->Barriers.me(EM, x))
plot!(range(0, C, 500), x->Barriers.me(AM, x))
plot!(0:Δ:C, [mec ; mec[end]], 
    linetype=:steppost, label="coarse, complete div.")


Y = fit(Histogram, xs, weights(EM.Ep), 0:Δ:C).weights
CM = Barriers.CoarseModel(X=X, Δ=Δs, s=sd, m=m, u=u, λ=0.0)

