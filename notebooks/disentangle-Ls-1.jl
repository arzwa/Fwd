
# L equally spaced loci with on a chromosome of a given map length with
# a fixed Ls. To what extent can we disentangle L and s? Increasing L
# and decreasing s, keeping Ls fixed, would not keep mₑ fixed, as the
# latter depens also on linkage.

using Random, Fwd, ProgressMeter, StatsBase, Serialization
using Plots; plotsdefault()
using Barriers

rng = Random.seed!(67)
Ls  = 0.2
u   = 1e-5
m   = 1e-3
C   = 0.5
NA  = 1
NB  = 1000
ngen = 5*10^5

LL = [5,10,20,30,40,50]
res = map(LL) do L
    s = Ls/L
    xs = collect((C/2L):(C/L):(C-C/2L))
    R   = LinearMap(C)
    AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
    nA = collect(1:NA)
    nB = collect(1:NB) .+ NA
    xA = [ ones(Int, L) for _=1:NA]
    xB = [zeros(Int, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, 
        gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, 
        gpm=MB, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    hmr = Fwd.hmrecrate(xs)
    mpop, ts = Fwd.simulate!(rng, mpop, init_ts(mpop), ngen)
    (hmr, s, mpop, ts)
end

models = map(LL) do L
    s = Ls/L
    xs = collect((C/2L):(C/L):(C-C/2L))
    R = Fwd.rec_matrix(xs)
    loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M)
    AM = AeschbacherModel(m, [-s for i=1:L], xs)
    EM, AM
end

ys = map(res) do X
    x, _, _, y = Fwd.diffdiv(X[5])
    x, y
end

xx = range(0, C, 500)
map(1:length(LL)) do i
    me_the = quadgk(x->Barriers.me(models[i][1], x), 0, C)[1]/C
    t_the  = 1/me_the + NA 
    plot(ys[i]..., color=:black, alpha=0.2, yscale=:log10, ylim=(NB,ngen))
    plot!(xx, x->1 / Barriers.me(models[i][1], x) + NA, color=:black)
    plot!(xx, x->1 / (m*Barriers.gff(models[i][2], x)) + NA, color=:salmon)
    hline!([t_the], color=:black)
end |> x->plot(x..., 
    layout=(length(LL),1), 
    size=(300,150*length(LL)),
    left_margin=9Plots.mm, ylabel="\$T_{AB}\$")


gs = map(1:length(LL)) do i
    x, y = ys[i]
    spans = [(x[i] - (i == 1 ? 0 : x[i-1])) for i=1:length(x)]
    tab_mean = sum(spans .* y)
    me_emp = 1 / (tab_mean - NA)
    me_the = quadgk(x->Barriers.me(models[i][1], x), 0, C)[1]/C
    L = LL[i]
    xs = collect((C/2L):(C/L):(C-C/2L))
    me_mid = Barriers.me(models[i][1], xs[L÷2])
    me_emp/m, me_the/m
end
plot(LL, exp.(-Ls ./ getindex.(res,2)), color=:black, marker=true, ms=3, label="\$g(\\overline{r})\$")
plot!(LL, getindex.(gs,2), color=:firebrick, marker=true, ms=3, label="\$\\overline{g}\$")
vline!(LL, color=:gray, alpha=0.2, ls=:dash, label="")
plot!(LL, first.(gs), ms=3, marker=true, ls=:dash, color=:gray,
    ylabel="\$g\$", xlabel="\$L\$", label="\$\\hat{g}\$")
plot!(legend=:topright, xlim=(0,LL[end]+5), size=(300,200))

getindex.(res,3) .* NB

plot(LL, getindex.(res,2))
plot!(LL, getindex.(res,3))
plot!(twinx(), LL, getindex.(res,2) ./ getindex.(res,3), 
    framestyle=:default, color=:black)

plot(LL, map(mean, last.(res)))


