@everywhere using Random, Fwd, ProgressMeter, StatsBase
using Plots, Barriers, Distributions; plotsdefault()

rng = Random.seed!(672)
L   = 10 
s   = 0.1/L
Ns  = 5.
u   = s/200
m   = s*1.0
r   = s/2
d   = Fwd.distance(r)
xs  = cumsum(fill(d, L))
C   = last(xs) + d
R   = LinearMap(C)
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
R   = LinearMap(C)
NA  = 1
NB  = ceil(Int, Ns/s)
nA  = collect(1:NA)
nB  = collect(1:NB) .+ NA
xA  = [ ones(Bool, L) for _=1:NA]
xB  = [zeros(Bool, L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
title = @sprintf "\$L=%d, Ls=%.2f, N_es=%.f, m/s=%.2f, r/s=%.2f\$" L L*s Ns m/s r/s

# Simulation
ngen = 2*10^5
nrep = 10
res  = pmap(1:nrep) do _
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
    every = ceil(Int, ngen/1000)
    mpop, ts, qs = simulate!(mpop, init_ts(mpop), ngen,
        x->mean(x.popB.x), every=every)
    seed, mpop, ts, qs
end

function tab(tss, nse=2)
    tbs = map(tss) do ts
        _ts = Fwd._add_grand_ancestor(ts)
        Fwd.diffdiv(_ts)[[1,4]]
    end
    xx, yy = Fwd.summarize_wins(tbs)
    ys = map(eachcol(yy)) do y
        # Assume coalescence times are Geometrically distributed with 
        # a noninformative Beta prior for the parameter of the
        # geometric distribution, determine the posterior Beta for the
        # parameter. Get [0.025, 0.975] posterior quantiles for the
        # Geometric distributions
        dp = Beta(length(y), sum(y) - 1)
        yu = mean(Geometric(quantile(dp, 0.025)))
        yl = mean(Geometric(quantile(dp, 0.975)))
        ym = 1/mean(dp)
        ym, yl, yu
    end
    ym = first.(ys)
    yl = getindex.(ys,2)
    yu = getindex.(ys,3)
    xx, ym, ym .- yl, yu .- ym 
end

xx, yb, yl, yu = tab(getindex.(res, 3))

Q = mapreduce(mean ∘ last, hcat, res)
pp = 1 .- mean(Q, dims=2) |> vec

P1 = plot(xx, yb, ribbon=(yl, yu), color=:lightgray, 
    label="", legend=:outertopright, 
    size=(700,300), title=title, 
    xlabel="map position (M)", 
    ylabel="\$T\$", margin=5Plots.mm)
BP = Equilibrium(BPModel(m, fill(s, L), xs, NB, u))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x),
    yscale=:log10, lw=1,label="BP, \$E[p]\$")
_BP = deepcopy(BP); _BP.Ep .= pp 
plot!(range(0, C, 500), x->1/Barriers.me(_BP, x),
    yscale=:log10, lw=1, label="BP, \$\\hat{p}\$")
loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
M = Equilibrium(Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1))
_BP = deepcopy(BP); _BP.Ep .= M.Ep 
plot!(range(0, C, 500), x->1/Barriers.me(_BP, x),
    yscale=:log10, lw=1, label="BP-ZSF24")
plot!(range(0, C, 500), x->1/Barriers.me(M, x), label="ZSF24", ls=:dot,
    yscale=:log10, lw=1)
AM = AeschbacherModel(m, fill(s, L), xs)
plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14",
    yscale=:log10, lw=1, ls=:dot)
hline!([ngen], ls=:dash, label="", color=:gray)

P2 = scatter(pp, M.Ep, label="ZSF24", legend=:topleft)
scatter!(pp, BP.Ep, label="BP")
plot!(x->x, ylabel="predicted", xlabel="simulation", 
    label="", xlim=(0,1), ylim=(0,1), color=:gray)

plot(P2, P1, layout=grid(1,2,widths=[0.3,0.7]), size=(800,270))


# Compare allele freq distributions
using WrightDistribution
QQ = hcat(vcat(map(X->X[2:end], last.(res))...)...)

k = 1
bpg = Barriers.gffs_(BP.model, BP.Ep)
mek = bpg[k]*BP.model.m
dk = Wright(-2Ns, NB*u, NB*(mek + u), 0.5)

bpg2 = Barriers.gffs_(BP.model, pp)
mek2 = bpg2[k]*BP.model.m
dk2 = Wright(-2Ns, NB*u, NB*(mek2 + u), 0.5)

stephist(QQ[k,:], xlabel="\$q\$", ylabel="density", legend=:topright, label="",
bins=0:0.02:1, norm=true, color=:black, fill=true, alpha=0.1, size=(300,200))
plot!(range(0,1,200), x->pdf(dk, 1-x), label="BP iteration") 
plot!(range(0,1,200), x->pdf(dk2, 1-x), label="BP empirical \$p\$") 


# mₑ point of view
plot(xx, 1 ./ yb, color=:lightgray)
plot!(range(0, C, 500), x->Barriers.me(BP, x), label="")
plot!(range(0, C, 500), x->Barriers.me(M, x), label="", ls=:dot)
plot!(range(0, C, 500), x->Barriers.me(AM, x), label="", ls=:dot)
_BP = deepcopy(BP); _BP.Ep .= pp 
plot!(range(0, C, 500), x->Barriers.me(_BP, x), label="")


# Check swamping threshold prediction
mss = range(0.05, 2.5, 100)
ps1 = map(mss) do ms
    BP = Equilibrium(BPModel(m=ms*s, s=fill(s, L), xs=xs, Ne=float(NB), u=u))
    BP.Ep
end |> x->hcat(x...)
ps2 = map(mss) do ms
    loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
    A = Barriers.Architecture(loci, xs, Fwd.rec_matrix(xs))
    M = Equilibrium(Barriers.MainlandIslandModel(arch=A, 
        m=ms*s, N=NB, mode=1))
    M.Ep
end |> x->hcat(x...)

i = 3
plot(mss, ps1[i,:])
plot!(mss, ps2[i,:])

ngen = 100_000
mss2 = range(extrema(mss)..., 10)
res2 = pmap(mss2) do ms
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    mpop = Fwd.TwoPopOneWay(ms*s, deepcopy(popA), deepcopy(popB))
    every = ceil(Int, ngen/1000)
    mpop, ts, qs = simulate!(mpop, init_ts(mpop), ngen,
        x->mean(x.popB.x), every=every)
    seed, mpop, ts, qs
end

P = 1 .- mapreduce(mean ∘ last, hcat, res2)
map([1,3,6,9]) do i
    plot(mss, ps1[i,:], title="locus $i", xlabel="\$m/s\$", ylabel="\$p\$", label="BP")
    plot!(mss, ps2[i,:], label="ZSF")
    scatter!(mss2, P[i,:], legend=i==1 ? :topright : false, 
        label="simulation", color=:black, ms=2)
    vline!([m/s], color=:lightgray, ls=:dash, label="")
end |> x->plot(x..., size=(500,400))


