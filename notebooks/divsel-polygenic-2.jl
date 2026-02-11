

using Distributed; addprocs(10)
@everywhere using Random, Fwd, ProgressMeter, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Barriers
using LinearAlgebra
import TreeSequences as TS

function estimate_coaltimes(tss, ci=0.95, pa=0., pb=0.; idx=4)
    q0 = (1-ci)/2
    q1 = 1-q0
    tbs = map(tss) do ts
       _ts = TS._add_grand_ancestor(ts)
       TS.diffdiv(_ts)[[1,idx]]
    end
    xx, yy = TS.summarize_wins(tbs)
    ys = map(eachcol(yy)) do y
       # Assume coalescence times are Geometrically distributed with
       # a noninformative Beta prior for the parameter of the
       # geometric distribution, determine the posterior Beta for the
       # parameter. Get [0.025, 0.975] posterior quantiles for the
       # Geometric distributions
       dp = Beta(length(y) + pa, sum(y) - 1 + pb)
       yu = mean(Geometric(quantile(dp, q0)))
       yl = mean(Geometric(quantile(dp, q1)))
       ym = 1/mean(dp)
       ym, yl, yu
    end
    ym = first.(ys)
    yl = getindex.(ys,2)
    yu = getindex.(ys,3)
    xx, ym, ym .- yl, yu .- ym
end


rng = Random.seed!(155)
Ls  = 0.1
L   = 25
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
u    = s̄/200
m    = s̄
NA   = 1
Ns   = 5. 
NB   = ceil(Int64, Ns/s̄)
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns m/s̄

map([0.1, 0.25, 0.5, 1.0]) do C
    xs = ys .* C
    R = Fwd.rec_matrix(xs)
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M);
    AM = AeschbacherModel(m, ss, xs)
    BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u))
    P = plot(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10)
    plot!(range(0, C, 500), x->1/Barriers.me(AM, x))
    plot!(range(0, C, 500), x->1/Barriers.me(BP, x))
    hline!([200000], ls=:dot, color=:lightgray)
end |> x->plot(x..., layout=(2,2), size=(800,400))
    
#ngen ~ 10 × max cross-pop coalescence time
#minme = quantile(map(x->Barriers.me(EM, x), range(0, C, 100)), 0.05)
#ngen = ceil(Int, (10 / minme * 1e-4)) * 10^4  
#ngen = 200_000
nrep = 10

res = map([0.25]) do C 
    ngen = 1_000_000
    res = pmap(1:nrep) do _
        xs  = ys .* C
        R   = LinearMap(C)
        AA  = Architecture([BiAllelic(0.0)   for _=1:L], xs, R)
        AB  = Architecture([BiAllelic(u)     for _=1:L], xs, R)
        MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
        MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
        nA  = collect(1:NA)
        nB  = collect(1:NB) .+ NA
        xA  = [ ones(Int, L) for _=1:NA]
        xB  = [zeros(Int, L) for _=1:NB]
        popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
        popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
        mpop = Fwd.TwoPopOneWay(m, popA, popB)
        pop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), ngen, 
            pop->mean(pop.popB.x), simplify=100, every=1000)
    end
end

serialize("data/2026-02-03-b.jls", res)

res = deserialize("data/2026-02-03-b.jls")
Ts = map(res) do X
    estimate_coaltimes(getindex.(X, 2), idx=4)
end

Qs = map(res) do X
    mean(mapreduce(mean, hcat, last.(X)), dims=2) |> vec
end

k = 1
Cs = [0.25,1.0]
plot(Ts[k][1], Ts[k][2], ribbon=Ts[k][3:4], yscale=:log10,
    framestyle=:default, color=:lightgray, ylabel="\$T\$", 
    margin=4Plots.mm, xlabel="map position (M)", xlim=(0,Cs[k]*1.01),)
#hline!([1/m], color=:lightgray, ls=:dash)
sticks!(twinx(), ys .* Cs[k], ss, framestyle=:default, ylim=(0,0.02),
    color=:firebrick, lw=3, alpha=0.3, ylabel="\$s\$", 
    xlabel="", title=title, size=(600,200), xlim=(0,Cs[k]*1.01))

Cs = [0.25,]
PP = map(enumerate(Ts)) do (i,(xx, yb, yl, yu))
    P1 = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8,
        label="", legend=:topright, 
        size=(600,500), title=i==1 ? title : "", 
        xlabel="map position (M)", yscale=:log10,
        ylabel="\$T\$", margin=5Plots.mm)
    C = Cs[i]
    xs = ys .* C
    R = Fwd.rec_matrix(xs)
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M);
    AM = AeschbacherModel(m, ss, xs)
    BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u))
    plot!(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10, alpha=0.8, label="ZSF24")
    plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14", alpha=0.8)
    plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP", alpha=0.8)
    hline!([200000], ls=:dot, color=:lightgray, label="")
    qs = Qs[i]
    P2 = scatter(qs, 1 .- EM.Ep, title="map length $(Cs[i])", ms=2, label="ZSF24")
    scatter!(qs, 1 .- BP.Ep, ms=2, label="BP", legend=:bottomright)
    plot!(x->x, color=:lightgray, ls=:dot, label="", 
        xlabel="\$\\hat{q}\$ (simulation)", 
        ylabel="\$\\mathbb{E}[q]\$", xlim=(0,1), ylim=(0,1))
    plot(P1, P2, layout=grid(1,2,widths=[0.75,0.25]))
end |> x->plot(x..., layout=(length(x),1), size=(700,200*length(x)))

Ps = let Ps = [plot() for i=1:L]
    map(1:1) do j
        C = Cs[j]
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u))
        g = Barriers.gff_sel(BP.model, BP.Ep)
        g2 = Barriers.gffs(EM)
        QQ = permutedims(hcat(map(x->hcat(x...)[:,2:end], last.(res[j]))...))
        col = j == 1 ? :black : :teal
        map(1:L) do k
            P = Ps[k]
            title!(@sprintf("\$s_{%d} = %.4f\$", k, ss[k]))
            stephist!(P, QQ[:,k], norm=true, bins=0:0.05:1.01, color=col,
                fill=true, alpha=0.2)
            d = Wright(2NB*ss[k], NB*(u + m*g[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->pdf(d, x), color=col) 
            d = Wright(2NB*ss[k], NB*(u + m*g2[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->pdf(d, x), color=col, ls=:dash)
        end
    end
    Ps
end
plot(Ps..., layout=(5,5), size=(900,700))

Ps = let Ps = [plot() for i=1:L]
    map(1:1) do j
        C = Cs[j]
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u))
        g = Barriers.gff_sel(BP.model, BP.Ep)
        g2 = Barriers.gffs(EM)
        QQ = permutedims(hcat(map(x->hcat(x...)[:,102:end], last.(res[j]))...))
        col = j == 1 ? :black : :teal
        map(1:L) do k
            P = Ps[k]
            title!(@sprintf("\$s_{%d} = %.4f\$", k, ss[k]))
            q = QQ[:,k]
            hist = normalize(fit(Histogram, q, 0:0.05:1.0))
            es = hist.edges[1]
            xx = [(es[i+1]+es[i])/2 for i=1:length(es)-1]
            scatter!(P, xx, log.(hist.weights), color=col, ms=1.5)
            vline!([mean(QQ[:,k])], lw=5, alpha=0.2, color=col)
            #stephist!(P, QQ[:,k], norm=true, bins=0:0.05:1.01, color=col,
            #    fill=true, alpha=0.2)
            d = Wright(2NB*ss[k], NB*(u + m*g[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->logpdf(d, x), color=col) 
            vline!([mean(d)], color=col)
            d = Wright(2NB*ss[k], NB*(u + m*g2[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->logpdf(d, x), color=col, ls=:dash)
            vline!([mean(d)], color=col, ls=:dash)
            y0, y1 = ylims(P)
            plot!(P, ylims=(-4, y1))
        end
    end
    Ps
end
P2 = plot(Ps..., layout=(5,5), size=(900,700))


Ps = let Ps = [plot() for i=1:L]
    map(1:1) do j
        C = Cs[j]
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u))
        g = Barriers.gff_sel(BP.model, BP.Ep)
        g2 = Barriers.gffs(EM)
        QQ = permutedims(hcat(map(x->hcat(x...)[:,102:end], last.(res[j]))...))
        col = j == 1 ? :black : :teal
        map(1:L) do k
            P = Ps[k]
            title!(@sprintf("\$s_{%d} = %.4f\$", k, ss[k]))
            q = QQ[:,k]
            hist = normalize(fit(Histogram, q, 0:0.05:1.0))
            es = hist.edges[1]
            xx = [(es[i+1]+es[i])/2 for i=1:length(es)-1]
            scatter!(P, xx, hist.weights, color=col, ms=2)
            vline!([mean(QQ[:,k])], lw=5, alpha=0.2, color=:black)
            #stephist!(P, QQ[:,k], norm=true, bins=0:0.05:1.01, color=col,
            #    fill=true, alpha=0.2)
            d = Wright(2NB*ss[k], NB*(u + m*g[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->pdf(d, x), color=col) 
            vline!([mean(d)], color=col)
            d = Wright(2NB*ss[k], NB*(u + m*g2[k]), NB*u, 0.5)
            plot!(P, range(0,1,200), x->pdf(d, x), color=col, ls=:dash)
            vline!([mean(d)], color=col, ls=:dash)
        end
    end
    Ps
end
plot(Ps..., layout=(5,5), size=(900,700))

# Coarse approximations ... (deprecated)
C  = 0.25 
nwin = 200
Δ  = step(range(0, C, nwin+1))
Δs = fill(Δ, nwin)
Xs = fit(Histogram, ys .* C, 0:Δ:C).weights
sw = fit(Histogram, ys .* C, weights(ss), 0:Δ:C).weights  
sw[isnan.(sw)] .= 0.0

CM = Barriers.CoarseModel2(X=Xs, Δ=Δs, m=m, u=u, s=sw, λ=1/NB) 
_gff, _CM = Barriers.gff(CM)
mec = _gff .* m

bs = 0:Δ:C-Δ
plot(bs, _gff, line=:steppost)
plot!(range(0, C, 500), x->Barriers.me(EM, x)/m, color=:black)
plot!(twinx(), bs, sw, line=:steppost, fill=true, alpha=0.5, color=:lightgray)

xy = map(X->Fwd.diffdiv(Fwd._add_grand_ancestor(X[2]), windows=collect(0:Δ:C))[[1,4]], res)
x, Y = Fwd.summarize_wins(first.(xy), last.(xy))
y = vec(mean(Y, dims=1))

plot(x, y, linetype=:steppre, color=:lightgray, alpha=1)
plot!(bs, 1 ./ mec, yscale=:log10, linetype=:steppost)
plot!(range(0, C, 500), x->1/Barriers.me(EM, x))

#CM = Barriers.CoarseModel(X=ones(length(Y)), Δ=Δs, s=Y, m=m, u=u, λ=0.)
#mec = m .* Barriers.gff(CM)
#plot!(bs, 1 ./ mec, linetype=:steppost, label="coarse 3")

zs = fit(Histogram, xs, weights(EM.Ep), 0:Δ:C).weights ./ Xs
Ep = _CM.Ep
Ep[Ep .== 0.0] .= NaN
plot(bs, _CM.Ep)
plot!(bs, zs)

sd = mean(ss .* EM.Ep)
plot(x, y, yscale=:log10, linetype=:steppre, color=:gray, alpha=0.5, label="")
plot!(range(0, C, 500), x->1/Barriers.me(EM, x), label="Zwaenepoel et al. 2024")
plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="Aeschbacher et al. 2017")
bs = 0:Δ:C-Δ
CM = Barriers.CoarseModel(X=Xs, Δ=Δs, s=sw, m=m, u=u, λ=0.)
mec = m .* Barriers.gff(CM)
plot!(bs, 1 ./ mec, linetype=:steppost, label="Approx. 1")
#Y = fit(Histogram, xs, weights(ss .* EM.Ep), 0:Δ:C).weights
#CM = Barriers.CoarseModel(X=ones(length(Y)), Δ=Δs, s=Y, m=m, u=u, λ=0.)
#mec = m .* Barriers.gff(CM)
#plot!(bs, 1 ./ mec, linetype=:steppost, label="hack")
CM = Barriers.CoarseModel2(X=Xs, Δ=Δs, m=m, u=u, s=sw, λ=1/NB) 
gff, _CM = Barriers.gff(CM)
mec = gff .* m
plot!(bs, 1 ./ mec, linetype=:steppost, label="Approx. 2")
plot!(legend=:outertopright, size=(1000,300))


# Y is Ls in windows
Y = fit(Histogram, xs, weights(ss), 0:Δ:C).weights
@assert sum(Y) ≈ sum(ss)

gg = map(1:nwin) do j
    loggj = 0.
    for i=1:nwin 
#        i == j && continue
        rij = Fwd.recrate(abs(i-j)*Δ)
        loggj -= Y[i]/(m + s̄ + rij)
    end
    loggj
end

P1 = plot(x, (1 ./ (y .- NA)) ./ m, linetype=:steppre, color=:gray, alpha=0.5)
plot!(Δ:Δ:C, exp.(gg), linetype=:steppre, color=:black)
plot!(bs, mec ./ m, linetype=:steppost)


# Zoom in


dd = map(X->Fwd.diffdiv(Fwd._add_grand_ancestor(X[2])), res)
xy = map(x->x[[1,4]], dd)
Q = mapreduce(x->vec(mean(x, dims=1)), hcat,getindex.(res, 3))
x, Y = Fwd.summarize_wins(first.(xy), last.(xy))
y = vec(mean(Y, dims=1))

CM = Barriers.CoarseModel2(X=Xs, Δ=Δs, m=m, u=u, s=sw, λ=1/NB) 
_gff, _CM = Barriers.gff(CM)
mec = _gff .* m
mm = mean(mec)

mesc(m, r, p) = m*r*(1-p)/(m*p + r*(1-p))

function mescx(x, m, Ep, xs)
    i = argmin(abs.(xs .- x))
    p = Ep[i]
    r = Fwd.recrate(abs(xs[i] - x))
    mesc(m, r, p)
end

function tc(x, m, Ep, xs; rmax=Inf)
    i = argmin(abs.(xs .- x))
    p = Ep[i]; q=1-p
    r = Fwd.recrate(abs(xs[i] - x))
    r <= rmax ? (1/r)*(p/q) + 1/m : 1/m
end


xspan = (0.0, 0.05) .+ 0.43 
plot(x,y, yscale=:log10, xlim=xspan, 
    framestyle=:default, color=:black)
plot!(bs, 1 ./ mec, color=:teal, line=:steppost)
hline!([1/mean(mec)], color=:teal, ls=:dash)
plot!(range(xspan..., 200), 
    x->tc(x, Barriers.me(EM,x), EM.Ep, xs), color=:orange,lw=2)
plot!(range(xspan..., 200), 
    x->1/Barriers.me(AM, x), color=:cyan,lw=2)
sticks!(twinx(), xs, EM.Ep, color=:firebrick, 
    alpha=0.5, lw=3, xlim=xspan, framestyle=:default)

xspan = (0.0, C)
plot(x,y, yscale=:log10, xlim=xspan, 
framestyle=:default, color=:black, size=(900,200))
plot!(bs, 1 ./ mec, color=:teal, line=:steppost)
hline!([1/mean(mec)], color=:teal, ls=:dash)
plot!(range(xspan..., 1000), 
    x->min(1e6, tc(x, Barriers.me(EM,x), EM.Ep, xs, rmax=s̄/4)), 
    color=:orange,lw=2,alpha=0.8)
plot!(range(xspan..., 1000), 
    x->1/Barriers.me(AM, x), color=:cyan,lw=1, ylim=(-Inf,1e6))
sticks!(twinx(), xs, EM.Ep, color=:firebrick, 
    alpha=0.1, lw=3, xlim=xspan, framestyle=:default)

xspan = (0.0, C)
xspan = (0.0, 0.05) .+ 0.43 
plot(x,y, yscale=:log10, xlim=xspan, ylabel="\$T\$", xlabel="map position (M)",
    framestyle=:default, color=:gray, size=(700,180), margin=5Plots.mm)
plot!(range(xspan..., 500), 
    x->1/Barriers.me(AM, x), color=1, alpha=0.7, lw=1.5, ylim=(-Inf,1e6))
plot!(range(xspan..., 500), 
    x->1/Barriers.me(BPE, x), color=2, alpha=0.7, lw=1.5, ylim=(-Inf,1e6))
plot!(range(xspan..., 500), 
    x->1/Barriers.me(EM, x), color=:black, ls=:dot, lw=2, ylim=(-Inf,1e6),)
sticks!(twinx(), xs, EM.Ep, color=:firebrick, ylabel="\$\\mathbb{E}[p]\$",
    alpha=0.2, lw=3, xlim=xspan, framestyle=:default, ylim=(0,1))


# simulation with decently sized mainland
rng = Random.seed!(25)
Ls  = 0.1
L   = 50
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
u   = s̄/200
m   = s̄/2
NA  = 1000
Ns  = 5. 
NB  = ceil(Int64, Ns/s̄)
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
C   = 0.5
G   = C*100*10^6   # (1cM/Mb)
xs  = ceil.(Int64, ys .* G)

BP = Equilibrium(BPModel(m=m, s=ss, xs=ys .* C, Ne=float(NB), u=u))
plot(range(0, C, 500), x->1/Barriers.me(BP, x), yscale=:log10)

R   = LinearPhysMap(C=C, G=G)
AA  = Architecture([BiAllelic(0.0)   for _=1:L], xs, R)
AB  = Architecture([BiAllelic(u)     for _=1:L], xs, R)
MA  = GPMap([HaploidLocus(0.0, i)    for i=1:L])
MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
nA  = collect(1:NA)
nB  = collect(1:NB) .+ NA
xA  = [ ones(Int, L) for _=1:NA]
xB  = [zeros(Int, L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 200_000
pop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), ngen, 
    pop->mean(pop.popB.x), simplify=100, every=1000)
    
BP = Equilibrium(BPModel(m=m, s=ss, xs=ys*C, Ne=float(NB), u=u))

xx, _,_,tt = TS.diffdiv(ts)
plot(xx, tt, color=:lightgray, yscale=:log10)
nn = 500 
tt = map(x->1/Barriers.me(BP, x), range(0, C, nn))
plot!(range(0, ts.L, nn), tt)





# -------------------
# Coarse model experiments. The data comes in windows, so what we need is
# an expected mₑ in a window (this is already a non-trivial approximation,
# substituting expected mₑ in a likelihood calculation instead of
# integrating the likelihood over mₑ or suchlike).
#
# CoarseModel2
#    X  :: Vector{Int}  # number of selected sites in window
#    Δ  :: Vector{T}    # winsizes in Morgan
#    R  :: Matrix{T} = winrecrates(Δ)  # between window recombination rates
#    s  :: Vector{T}  # selection coefficient/vector of selection coefficients
#    Ep :: Vector{T} = ones(length(X))
#    m  :: T  # migration rate
#    u  :: T  # mutation rate
#    λ  :: T  # coal. rate (inverse pop size)


C  = 0.25 
nwin = 200
ws = range(0, C, nwin+1)
Δ  = step(ws)
Δs = fill(Δ, nwin)
xs = ys .* C
Xs = fit(Histogram, xs, 0:Δ:C).weights
sw = fit(Histogram, xs, weights(ss), 0:Δ:C).weights  ./ Xs
sw[isnan.(sw)] .= 0.0

CM = Barriers.CoarseModel2(X=Xs, Δ=Δs, s=sw, m=m, u=u, λ=1/NB)
gs, CM_ = Barriers.gff(CM)
plot(ws[2:end], 1 ./ (m .* gs), yscale=:log10)
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
EM = Barriers.Equilibrium(M);
AM = AeschbacherModel(m, ss, xs)
BP = Equilibrium(BPModel(m, ss, xs, NB, u))
plot(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10, label="ZSF24")
plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14")
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")

BP = Equilibrium(BPModel(m, ss, xs, NB, u))
plot(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP", yscale=:log10)
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=1))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=2))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=4))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=8))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=16))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
BP = Equilibrium(BPModel(m, ss, xs, NB, u, n=32))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")


# Swamping etc.
rng = Random.seed!(155)
Ls  = 0.1
L   = 25
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
C   = 0.10
xs  = ys .* C
u    = s̄/200
NA   = 1
Ns   = 5. 
NB   = ceil(Int64, Ns/s̄)
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns m/s̄

mss = range(0.05, 2.0, 40)
pps = map(mss) do ms
    m = ms*s̄
    R = Fwd.rec_matrix(xs)
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M);
    BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u))
    EM.Ep, BP.Ep
end
P1 = hcat(first.(pps)...)
P2 = hcat(last.(pps)...)

map(1:L) do k
    plot(mss, P1[k,:])
    plot!(mss, P2[k,:])
end |> x->plot(x..., layout=(5, L÷5), size=(160L÷5, 5*140))


