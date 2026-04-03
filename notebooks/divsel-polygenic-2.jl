
using Distributed#; addprocs(10)
@everywhere using Random, Fwd, ProgressMeter, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Barriers
using LinearAlgebra
using Parameters
import TreeSequences as TS
import MCMCChains: mcse

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


rng = Random.seed!(135)
Ls  = 0.25
L   = 50
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
u    = s̄/200
NA   = 1
Ns   = 5. 
NB   = ceil(Int64, Ns/s̄)
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]

# C = 0.25 is very tight linkage...
C = 0.5
Fwd.rbar(ys * C)

mss = [1, 1.5, 2]
map(mss) do ms
    m = ms*s̄
    map([0.25, 0.5, 1.0]) do C
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        AM = AeschbacherModel(m, ss, xs)
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u), α=0.2)
        P = plot(range(0, C, 500),  x->1/Barriers.me(EM, x), title="\$C=$C, m/\\bar{s}=$ms\$", 
            label="ZSF24", yscale=:log10, legend=C==1.0 ? :topright : false)
        plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14")
        plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP")
        plot!(ylim=(50,3*10^6), xlabel="map position", 
            ylabel="\$T\$", left_margin=3Plots.mm, bottom_margin=5Plots.mm)
    #    hline!([200000], ls=:dot, color=:lightgray, label="", ylim=(1000, 2*10^6))
    end |> x->plot(x..., layout=(1,3), size=(900,220))
end |> x->plot(x..., layout=(length(mss), 1), size=(900,220*length(mss)))

Fwd.rec_matrix(ys .* 0.25) ./ s̄
    
#ngen ~ 10 × max cross-pop coalescence time
#minme = quantile(map(x->Barriers.me(EM, x), range(0, C, 100)), 0.05)
#ngen = ceil(Int, (10 / minme * 1e-4)) * 10^4  
#ngen = 200_000
nrep = 5
ms = 1.5
m = ms*s̄
res = map([0.25, 0.5, 1.0]) do C 
    ngen = 2_000_000
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

#serialize("data/2026-03-05-ms1.5.jls", res)

data = [
    (ms=1.0, res=deserialize("data/2026-03-05-ms1.0.jls")),
    (ms=1.5, res=deserialize("data/2026-03-05-ms1.5.jls"))
]

Tss = map(data) do x
    Ts = map(x.res) do X
        estimate_coaltimes(getindex.(X, 2), idx=4)
    end
end

Tsw = map(data) do x
    Ts = map(x.res) do X
        estimate_coaltimes(getindex.(X, 2), idx=3)
    end
end

Qss = map(data) do x
    Qs = map(x.res) do X
        Q = map(last.(X)) do ys
            Q = hcat(ys...)
        end |> Q->hcat(Q...)
        se = mcse.(eachrow(Q))
        mn = vec(mean(Q, dims=2))
        mn, se
    end
end

k = 1
ms = data[k].ms
m = ms*s̄
Cs = [0.25,0.5,1]
title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns ms
PP = map(enumerate(Tss[k])) do (i,(xx, yb, yl, yu))
    P1 = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8,
        label="", legend=i==1 ? :topleft : false, 
        size=(600,500), title=i==1 ? title : "", 
        xlabel="map position (M)", yscale=:log10,
        ylabel="\$T\$", margin=1Plots.mm)
    C = Cs[i]
    xs = ys .* C
    R = Fwd.rec_matrix(xs)
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M);
    AM = AeschbacherModel(m, ss, xs)
    BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u), α=0.2)
    plot!(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10, alpha=0.8, label="ZSF24")
    plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14", alpha=0.8)
    plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP", alpha=0.8)
    qs, se = Qss[1][i]
    P2 = scatter(qs, 1 .- EM.Ep, title="map length $(Cs[i])", ms=2, label="ZSF24")
    scatter!(qs, 1 .- BP.Ep, ms=2, label="BP", legend=:bottomright)
#    P2 = scatter(qs, 1 .- EM.Ep, xerr=3se, msc=1, lw=1, title="map length $(Cs[i])", ms=2, label="ZSF24")
#    scatter!(qs, 1 .- BP.Ep, ms=2, xerr=3se, msc=2, lw=1, label="BP", legend=:bottomright)
    plot!(x->x, color=:lightgray, ls=:dot, label="", 
        xlabel="\$\\hat{q}\$ (simulation)", 
        ylabel="\$\\mathbb{E}[q]\$", xlim=(0,1), ylim=(0,1))
    plot(P1, P2, layout=grid(1,2,widths=[0.8,0.2]))
end |> x->plot(x..., layout=(length(x),1), size=(750,200*length(x)))

# within pop coal time
k = 1
ms = data[k].ms
m = ms*s̄
Cs = [0.25,0.5,1]
title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns ms
PP = map(enumerate(Tsw[k])) do (i,(xx, yb, yl, yu))
    P1 = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8,
        label="", legend=i==1 ? :topleft : false, 
        size=(600,500), title=i==1 ? title : "", 
        xlabel="map position (M)", yscale=:log10,
        ylabel="\$T\$", margin=1Plots.mm)
    C = Cs[i]
    xs = ys .* C
    R = Fwd.rec_matrix(xs)
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M);
    AM = AeschbacherModel(m, ss, xs)
    BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u), α=0.2)
    tw(me) = NB*(3 + 2me*NA)/(1 + 2me*NB)
    plot!(range(0, C, 500),  x->tw(Barriers.me(EM, x)), yscale=:log10, alpha=0.8, label="ZSF24")
    plot!(range(0, C, 500), x->tw(Barriers.me(AM, x)), label="AB14", alpha=0.8)
    plot!(range(0, C, 500), x->tw(Barriers.me(BP, x)), label="BP", alpha=0.8)
    qs, se = Qss[1][i]
    P2 = scatter(qs, 1 .- EM.Ep, title="map length $(Cs[i])", ms=2, label="ZSF24")
    scatter!(qs, 1 .- BP.Ep, ms=2, label="BP", legend=:bottomright)
#    P2 = scatter(qs, 1 .- EM.Ep, xerr=3se, msc=1, lw=1, title="map length $(Cs[i])", ms=2, label="ZSF24")
#    scatter!(qs, 1 .- BP.Ep, ms=2, xerr=3se, msc=2, lw=1, label="BP", legend=:bottomright)
    plot!(x->x, color=:lightgray, ls=:dot, label="", 
        xlabel="\$\\hat{q}\$ (simulation)", 
        ylabel="\$\\mathbb{E}[q]\$", xlim=(0,1), ylim=(0,1))
    plot(P1, P2, layout=grid(1,2,widths=[0.8,0.2]))
end |> x->plot(x..., layout=(length(x),1), size=(750,200*length(x)))

# with focus on a region
PP = map(enumerate(data)) do (j,Y)
    Ts = Tss[j]
    Qs = Qss[j]
    ms = Y.ms
    m = ms*s̄
    Cs = [0.25,0.5,1]
    title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns ms
    map(enumerate(Ts)) do (i,(xx, yb, yl, yu))
        P1 = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8, lw=0, la=0.,
            label="", legend=i==1 ? :topleft : false, 
            size=(600,500), title=i==1 ? title : "", 
            xlabel="map position (M)", yscale=:log10,
            ylabel="\$T\$", margin=1Plots.mm)
        C = Cs[i]
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        AM = AeschbacherModel(m, ss, xs)
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u), α=0.2)
        plot!(range(0, C, 500),  x->1/Barriers.me(EM, x), yscale=:log10, alpha=0.8, label="ZSF24")
        plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="AB14", alpha=0.8)
        plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="BP", alpha=0.8)
        qs, se = Qs[i]
        x0, x1 = (0.105, 0.19) .* C
        P1b = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8, lw=0, la=0.,
            xlabel="map position (M)", yscale=:log10, xlim=(x0,x1))
        plot!(range(x0, x1, 200),  x->1/Barriers.me(EM, x), yscale=:log10, alpha=0.8, label="ZSF24")
        plot!(range(x0, x1, 200), x->1/Barriers.me(AM, x), label="AB14", alpha=0.8)
        plot!(range(x0, x1, 200), x->1/Barriers.me(BP, x), label="BP", alpha=0.8)
        P2 = scatter(qs, 1 .- EM.Ep, title="map length $(Cs[i])", ms=2, label="ZSF24")
        scatter!(qs, 1 .- BP.Ep, ms=2, label="BP", legend=i == 1 ? :bottomright : false)
        plot!(x->x, color=:lightgray, ls=:dot, label="", 
            xlabel="\$\\hat{q}\$ (simulation)", 
            ylabel="\$\\mathbb{E}[q]\$", xlim=(0,1), ylim=(0,1))
        plot(P1, P1b, P2, layout=grid(1,3,widths=[0.6,0.2,0.2]))
    end |> x->plot(x..., layout=(length(x),1), size=(850,200*length(x)))
end 

PP = map(enumerate(data)) do (j,Y)
    Ts = Tsw[j]
    Qs = Qss[j]
    ms = Y.ms
    m = ms*s̄
    Cs = [0.25,0.5,1]
    title = @sprintf "\$L=%d, L\\bar{s}=%.2f, N_e\\bar{s}=%.f, m/\\bar{s}=%.2f\$" L Ls Ns ms
    map(enumerate(Ts)) do (i,(xx, yb, yl, yu))
        P1 = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8, lw=0, la=0.,
            label="", legend=i==1 ? :topleft : false, 
            size=(600,500), title=i==1 ? title : "", 
            xlabel="map position (M)", yscale=:log10,
            ylabel="\$T\$", margin=1Plots.mm)
        C = Cs[i]
        xs = ys .* C
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M);
        AM = AeschbacherModel(m, ss, xs)
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=NB, u=u), α=0.2)
        tw(me) = NB*(3 + 2me*NA)/(1 + 2me*NB)
        plot!(range(0, C, 500), x->tw(Barriers.me(EM, x)), yscale=:log10, alpha=0.8, label="ZSF24")
        plot!(range(0, C, 500), x->tw(Barriers.me(AM, x)), label="AB14", alpha=0.8)
        plot!(range(0, C, 500), x->tw(Barriers.me(BP, x)), label="BP", alpha=0.8)
        qs, se = Qs[i]
        x0, x1 = (0.105, 0.19) .* C
        P1b = plot(xx, yb, ribbon=(yl, yu), color=:gray, alpha=0.8, lw=0, la=0.,
            xlabel="map position (M)", yscale=:log10, xlim=(x0,x1))
        plot!(range(x0, x1, 200),  x->tw(Barriers.me(EM, x)), yscale=:log10, alpha=0.8, label="ZSF24")
        plot!(range(x0, x1, 200), x->tw(Barriers.me(AM, x)), label="AB14", alpha=0.8)
        plot!(range(x0, x1, 200), x->tw(Barriers.me(BP, x)), label="BP", alpha=0.8)
        P2 = scatter(qs, 1 .- EM.Ep, title="map length $(Cs[i])", ms=2, label="ZSF24")
        scatter!(qs, 1 .- BP.Ep, ms=2, label="BP", legend=i == 1 ? :bottomright : false)
        plot!(x->x, color=:lightgray, ls=:dot, label="", 
            xlabel="\$\\hat{q}\$ (simulation)", 
            ylabel="\$\\mathbb{E}[q]\$", xlim=(0,1), ylim=(0,1))
        plot(P1, P1b, P2, layout=grid(1,3,widths=[0.6,0.2,0.2]))
    end |> x->plot(x..., layout=(length(x),1), size=(850,200*length(x)))
end 

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
            plot!(P, ylims=(-8, y1))
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

tb = map(range(0, C, 500)) do x
    me = Barriers.me(BP, x)
    tab = 1/me + NA
    tb = 1/(2me + 1/NB) + 2me/(2me + 1/NB)*tab
    tb_ = NB*(3 + 2me*NA)/(1 + 2me*NB)
    @assert tb ≈ tb_
    x, tab, tb
end

plot(first.(tb), getindex.(tb,2), yscale=:log10)
plot!(first.(tb), getindex.(tb,3))

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


# Within pop coalescence times
# ==========================================================================
rng = Random.seed!(135)
Ls  = 0.25
L   = 50
s̄   = Ls/L
dfe = Exponential(s̄)
ss  = rand(rng, dfe, L)
ss .*= s̄/mean(ss) 
α   = 2.0
zs  = [0.0 ; cumsum(rand(rng, Dirichlet(L, α)))] 
ys  = [(zs[i] + zs[i+1])/2 for i=1:L]
u    = s̄/200
NA   = 1
Ns   = 5. 
NB   = ceil(Int64, Ns/s̄)
ms   = 1.0
m    = ms*s̄
C    = 0.25
xs   = ys .* C
BP   = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u), α=0.2)
xx, yy, y1, y2 = Tsw[1][1]

plot(xx, yy, ribbon=(y1,y2), color=:lightgray, size=(700,200), yscale=:log10)

fs, ws = Barriers.bc_fitnesses(BP, 50)
kmax = findfirst(k->1 - ws[k] < 1/NB, 1:length(ws))
chi = Barriers.bc_ancestries(BP, kmax) |> sum

Barriers.bc_ancestries(BP, 50) |> cumsum |> plot
plot!(twinx(), ws, color=2)

_tab(m) = 1/m + NA
_tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
pred = map(range(0, C, 500)) do x
    me = Barriers.me(BP, x)
    tb0 = _tw(NA, NB, me)
    tb1 = (1-2chi)*_tw(NA,NB,me) + 2chi*_tab(me) 
    x, tb0, tb1
end
plot(xx, yy, ribbon=(y1,y2), color=:lightgray, size=(700,200), label="", yscale=:log10)
plot!(getindex.(pred, Ref([1,2])), label="\$T_B(m_e)\$")
plot!(getindex.(pred, Ref([1,3])), label="veller")



# BP approach
# -----------
function chi2(m, r1, r2, sp1, sp2)
    m*(-r1*r2 + (sp1 + r1 + r2)*(sp2 + r1 + r2))/((sp1 + r1)*(sp2 + r2)*(sp1 + sp2 + r1 + r2))
end

_tab(m) = 1/m + NA
_tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
pred = map(range(0, C, 500)) do x
    me = Barriers.me(BP, x)
    r = Fwd.recrate(minimum(abs.(BP.model.xs .- x)))
    xx1 = me/r
    li = findlast(i->xs[i] <= x, 1:length(xs))
    xx = if isnothing(li)
        r = Fwd.recrate(xs[1] - x)
        me/r
    elseif li == length(xs)
        r = Fwd.recrate(x - xs[end])
        me/r
    else
        ri = li + 1
        rl = Fwd.recrate(x - xs[li])
        rr = Fwd.recrate(xs[ri] - x)
        xx = chi2(m, rl, rr, ss[li]*BP.Ep[li], ss[ri]*BP.Ep[ri])
        xx
        me*(rr + rl)/(rr*rl)
        #me/(√(rr*rl))
    end
    x, xx1, xx, (1-2xx)*_tw(NA,NB,me) + 2xx*_tab(me) , _tw(NA,NB,me)
end

plot(xx, yy, ribbon=(y1,y2), color=:lightgray, size=(700,200), yscale=:log10)
plot!(getindex.(pred, Ref([1,4])))
plot!(getindex.(pred, Ref([1,5])))

function rho(BP, kmax=floor(Int, log2(0.5/BP.model.m)))
    # heuristic for kmax: at most 50% non-residents
    @unpack m, s, xs = BP.model
    ws = map(0:kmax-1) do k
        exp(-sum([s[i]*BP.Ep[i] for i=1:length(xs)])/2^k)
    end 
#    fs = m*cumprod(ws)
    # fs = [m*W₀, m*W₀*W₁, ...]
    # i.e. [F1s,  BC1    , ...]
    fs = m*cumprod(2 .* ws)
#    return fs, ws
    # fs = [2m*W₀, 4m*W₀*W₁, ...]
    # i.e. [F1s,  BC1    , ...]
    # when we sample, there are never pure migrants around: we sample after
    # a migration + reproduction cycle. All migrants mate with residents to
    # make F1s.
    # The total probability to sample a lineage that traces back to a
    # migrant in the recent past ≈ (1/2)*fs[1] + (1/4)*fs[2] + ...
    sum([1/2^k * fs[k] for k=1:kmax])
    # If migrants were included and
    # fs = [M, F1, BC1, ...]
    # then we'd have something like: 
    # sum([1/2^k * fs[k+1] for k=0:kmax-1])
    # I guess we could make various corrections to our estimates for the fse
end

using Interpolations
xs_ = [0 ; xs ; C]
ys_ = [BP.Ep[1]; BP.Ep; BP.Ep[end]]
itp = linear_interpolation(xs_, ys_)
plot(range(0,C,500), x->itp(x))

tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
ρ = rho(BP, 5)
tb = map(range(0,C,500)) do x
    m = BP.model.m
    me = Barriers.me(BP, x)
    g = me/m
    q = 1-itp(x)
    tab = 1/me + NA
    tb = tw(NA, NB, me)
    x, (1-ρ)^2*tb + 2ρ*(1-ρ)*tab + ρ^2*NA, (1-2q)*tb + 2q*tab, tb
end
plot(xx, yy, ribbon=(y1,y2), color=:lightgray, size=(700,200))
plot!(getindex.(tb, Ref([1,2])), yscale=:log10, lw=2)
#plot!(getindex.(tb, Ref([1,3])), yscale=:log10, lw=2)
plot!(getindex.(tb, Ref([1,4])), yscale=:log10, lw=2)

plot!(first.(tb), x->1/Barriers.me(BP,x))

# Fit a model
using Optim
function objective(xs, ts, BP, NA, NB, trans=log)
    mes = map(x->Barriers.me(BP, x), xs)
    tab = 1 ./ mes .+ NA
    tb  = NB .* (3 .+ mes .* 2NA) ./ (1 .+ mes .* 2NB)
    return function obj(α)
        tpred = tab*α + (1-α)*tb
        tpred, sum((trans.(ts) .- trans.(tpred)) .^ 2)
    end
end

ps = map(1:2) do j
    ms = [1, 1.5][j]
    map(1:3) do k
        C = [0.25,0.5,1.0][k]
        xx, yy, y1, y2 = Tsw[j][k]
        m = ms*s̄
        xs = ys .* C
        BP = Equilibrium(BPModel(m=m, s=ss, xs=xs, Ne=float(NB), u=u), α=0.2)
        itp = linear_interpolation([0; xx; C], [yy[1] ; yy; yy[end]]) 
        n = 500
        xk = range(0,C,n)
        tt = [itp(x) for x in xk]
        f1 = objective(range(0, C, 500), tt, BP, NA, NB, x->1/x)
        result = Optim.optimize(x->f1(x)[2], 0, 1)
        α1 = result.minimizer
        plot(xx, yy, ribbon=(y1,y2), color=:lightgray, size=(700,200),
            label="")
#        plot!(xk, tt)
        g = mean(1 ./ (f1.tab .- NA)) /m
        ρ = rho(BP)
        plot!(xk, f1(α1)[1], yscale=:log10, 
            label=@sprintf("\$\\alpha = %.3f\$", α1),
            title=@sprintf("\$C=%.2f, m/s=%.2f\$", C, ms))
        tb = map(range(0,C,500)) do x
            me = Barriers.me(BP, x)
            tab = 1/me + NA
            tb = tw(NA, NB, me)
            x, (1-ρ)^2*tb + 2ρ*(1-ρ)*tab + ρ^2*NA, tb
        end
        plot!(getindex.(tb, Ref([1,2])), label=@sprintf("\$\\rho=%.3f\$", ρ))
        plot!(getindex.(tb, Ref([1,3])), label="\$m_e\$")
    end |> x->plot(x..., layout=(3,1), legend=:topright)
end 
plot(ps..., layout=(1,2), size=(900,600))

