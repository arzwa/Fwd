using Distributed
@everywhere using Fwd, TreeSequences, StatsBase, WrightDistribution, Barriers, TwoLocusModels
using Plots; plotsdefault()
using Serialization

@everywhere function twolocus_cb(pop, i, j; states=[[0,0],[0,1],[1,0],[1,1]])
    pm = proportionmap(map(x->x[[i,j]], pop.x))
    [haskey(pm, x) ? pm[x] : 0.0 for x in states]
end

function scpred(m, xs, D, q_1, q_2, x)
    if x < xs[1] || x > xs[2]
        r_1 = Fwd.recrate(min(abs(xs[1] - x), abs(xs[2]-x)))
        r = Fwd.recrate(xs[2] - xs[1])
        (-D^2*m*r_1 - D^2*r*r_1 + D^2*r_1^2 - D*m*q_1*q_2*r - 2*D*m*q_1*q_2*r_1 + D*m*q_1*r + D*m*q_1*r_1 + D*m*q_2*r - D*m*r + D*m*r_1 + 2*D*q_1*q_2*r_1^2 - D*q_1*r_1^2 - m*q_1^2*q_2^2*r - m*q_1^2*q_2^2*r_1 + m*q_1^2*q_2*r + m*q_1^2*q_2*r_1 + m*q_1*q_2^2*r - m*q_1*q_2*r + m*q_1*q_2*r_1 - m*q_1*r_1 + q_1^2*q_2^2*r*r_1 + q_1^2*q_2^2*r_1^2 - q_1^2*q_2*r*r_1 - q_1^2*q_2*r_1^2)/(m*r_1*(D^2*r_1 + D*q_1*q_2*r + 2*D*q_1*q_2*r_1 - D*q_1*r - D*q_1*r_1 + q_1^2*q_2^2*r + q_1^2*q_2^2*r_1 - q_1^2*q_2*r - q_1^2*q_2*r_1))
    else
        # in between
        r_1 = Fwd.recrate(x-xs[1])
        r_2 = Fwd.recrate(xs[2]-x)
        (D^4*m*r_1*r_2 + D^4*r_1^2*r_2 + D^4*r_1*r_2^2 - D^3*m*q_1*q_2*r_1^2 - D^3*m*q_1*q_2*r_2^2 + D^3*m*q_1*r_1^2 + D^3*m*q_1*r_1*r_2 + D^3*m*q_1*r_2^2 + D^3*m*q_2*r_1^2 + D^3*m*q_2*r_1*r_2 + D^3*m*q_2*r_2^2 - D^3*m*r_1^2 - 2*D^3*m*r_1*r_2 - D^3*m*r_2^2 - D^2*m*q_1^2*q_2^2*r_1^2 - 2*D^2*m*q_1^2*q_2^2*r_1*r_2 - D^2*m*q_1^2*q_2^2*r_2^2 + D^2*m*q_1^2*q_2*r_1^2 + 3*D^2*m*q_1^2*q_2*r_1*r_2 + D^2*m*q_1^2*q_2*r_2^2 - D^2*m*q_1^2*r_1*r_2 + D^2*m*q_1*q_2^2*r_1^2 + 3*D^2*m*q_1*q_2^2*r_1*r_2 + D^2*m*q_1*q_2^2*r_2^2 - D^2*m*q_1*q_2*r_1^2 - 5*D^2*m*q_1*q_2*r_1*r_2 - D^2*m*q_1*q_2*r_2^2 + 2*D^2*m*q_1*r_1*r_2 - D^2*m*q_2^2*r_1*r_2 + 2*D^2*m*q_2*r_1*r_2 - D^2*m*r_1*r_2 + D*m*q_1^3*q_2^3*r_1^2 + D*m*q_1^3*q_2^3*r_2^2 - 3*D*m*q_1^3*q_2^2*r_1^2 - 2*D*m*q_1^3*q_2^2*r_2^2 + 3*D*m*q_1^3*q_2*r_1^2 + D*m*q_1^3*q_2*r_2^2 - D*m*q_1^3*r_1^2 - 2*D*m*q_1^2*q_2^3*r_1^2 - 3*D*m*q_1^2*q_2^3*r_2^2 + 6*D*m*q_1^2*q_2^2*r_1^2 + 6*D*m*q_1^2*q_2^2*r_2^2 - 6*D*m*q_1^2*q_2*r_1^2 - 3*D*m*q_1^2*q_2*r_2^2 + 2*D*m*q_1^2*r_1^2 + D*m*q_1*q_2^3*r_1^2 + 3*D*m*q_1*q_2^3*r_2^2 - 3*D*m*q_1*q_2^2*r_1^2 - 6*D*m*q_1*q_2^2*r_2^2 + 3*D*m*q_1*q_2*r_1^2 + 3*D*m*q_1*q_2*r_2^2 - D*m*q_1*r_1^2 - D*m*q_2^3*r_2^2 + 2*D*m*q_2^2*r_2^2 - D*m*q_2*r_2^2 + m*q_1^4*q_2^4*r_1^2 + m*q_1^4*q_2^4*r_1*r_2 + m*q_1^4*q_2^4*r_2^2 - 3*m*q_1^4*q_2^3*r_1^2 - 2*m*q_1^4*q_2^3*r_1*r_2 - 2*m*q_1^4*q_2^3*r_2^2 + 3*m*q_1^4*q_2^2*r_1^2 + m*q_1^4*q_2^2*r_1*r_2 + m*q_1^4*q_2^2*r_2^2 - m*q_1^4*q_2*r_1^2 - 2*m*q_1^3*q_2^4*r_1^2 - 2*m*q_1^3*q_2^4*r_1*r_2 - 3*m*q_1^3*q_2^4*r_2^2 + 6*m*q_1^3*q_2^3*r_1^2 + 3*m*q_1^3*q_2^3*r_1*r_2 + 6*m*q_1^3*q_2^3*r_2^2 - 6*m*q_1^3*q_2^2*r_1^2 - 3*m*q_1^3*q_2^2*r_2^2 + 2*m*q_1^3*q_2*r_1^2 - m*q_1^3*q_2*r_1*r_2 + m*q_1^2*q_2^4*r_1^2 + m*q_1^2*q_2^4*r_1*r_2 + 3*m*q_1^2*q_2^4*r_2^2 - 3*m*q_1^2*q_2^3*r_1^2 - 6*m*q_1^2*q_2^3*r_2^2 + 3*m*q_1^2*q_2^2*r_1^2 - 3*m*q_1^2*q_2^2*r_1*r_2 + 3*m*q_1^2*q_2^2*r_2^2 - m*q_1^2*q_2*r_1^2 + 2*m*q_1^2*q_2*r_1*r_2 - m*q_1*q_2^4*r_2^2 - m*q_1*q_2^3*r_1*r_2 + 2*m*q_1*q_2^3*r_2^2 + 2*m*q_1*q_2^2*r_1*r_2 - m*q_1*q_2^2*r_2^2 - m*q_1*q_2*r_1*r_2 - q_1^4*q_2^4*r_1^2*r_2 - q_1^4*q_2^4*r_1*r_2^2 + 2*q_1^4*q_2^3*r_1^2*r_2 + 2*q_1^4*q_2^3*r_1*r_2^2 - q_1^4*q_2^2*r_1^2*r_2 - q_1^4*q_2^2*r_1*r_2^2 + 2*q_1^3*q_2^4*r_1^2*r_2 + 2*q_1^3*q_2^4*r_1*r_2^2 - 4*q_1^3*q_2^3*r_1^2*r_2 - 4*q_1^3*q_2^3*r_1*r_2^2 + 2*q_1^3*q_2^2*r_1^2*r_2 + 2*q_1^3*q_2^2*r_1*r_2^2 - q_1^2*q_2^4*r_1^2*r_2 - q_1^2*q_2^4*r_1*r_2^2 + 2*q_1^2*q_2^3*r_1^2*r_2 + 2*q_1^2*q_2^3*r_1*r_2^2 - q_1^2*q_2^2*r_1^2*r_2 - q_1^2*q_2^2*r_1*r_2^2)/(m*r_1*r_2*(D^3*q_1*q_2*r_1 + D^3*q_1*q_2*r_2 - D^3*q_1*r_1 - D^3*q_1*r_2 - D^3*q_2*r_1 - D^3*q_2*r_2 + D^3*r_1 + D^3*r_2 + D^2*q_1^2*q_2^2*r_1 + D^2*q_1^2*q_2^2*r_2 - D^2*q_1^2*q_2*r_1 - D^2*q_1^2*q_2*r_2 - D^2*q_1*q_2^2*r_1 - D^2*q_1*q_2^2*r_2 + D^2*q_1*q_2*r_1 + D^2*q_1*q_2*r_2 - D*q_1^3*q_2^3*r_1 - D*q_1^3*q_2^3*r_2 + 2*D*q_1^3*q_2^2*r_1 + 2*D*q_1^3*q_2^2*r_2 - D*q_1^3*q_2*r_1 - D*q_1^3*q_2*r_2 + 2*D*q_1^2*q_2^3*r_1 + 2*D*q_1^2*q_2^3*r_2 - 4*D*q_1^2*q_2^2*r_1 - 4*D*q_1^2*q_2^2*r_2 + 2*D*q_1^2*q_2*r_1 + 2*D*q_1^2*q_2*r_2 - D*q_1*q_2^3*r_1 - D*q_1*q_2^3*r_2 + 2*D*q_1*q_2^2*r_1 + 2*D*q_1*q_2^2*r_2 - D*q_1*q_2*r_1 - D*q_1*q_2*r_2 - q_1^4*q_2^4*r_1 - q_1^4*q_2^4*r_2 + 2*q_1^4*q_2^3*r_1 + 2*q_1^4*q_2^3*r_2 - q_1^4*q_2^2*r_1 - q_1^4*q_2^2*r_2 + 2*q_1^3*q_2^4*r_1 + 2*q_1^3*q_2^4*r_2 - 4*q_1^3*q_2^3*r_1 - 4*q_1^3*q_2^3*r_2 + 2*q_1^3*q_2^2*r_1 + 2*q_1^3*q_2^2*r_2 - q_1^2*q_2^4*r_1 - q_1^2*q_2^4*r_2 + 2*q_1^2*q_2^3*r_1 + 2*q_1^2*q_2^3*r_2 - q_1^2*q_2^2*r_1 - q_1^2*q_2^2*r_2))
    end
end

function scqle(m, xs, q_1, q_2, s1, s2, x)
    r = Fwd.recrate(xs[2]-xs[1])
    D = m*(1-q_1)*(1-q_2)/(m - 2*q_1*s1 - 2*q_2*s2 + r + s1 + s2)
    scpred(m, xs, D, q_1, q_2, x)
end

# Exploration
s   = 0.02
L   = 2
NB  = 500
s   = 10/NB
u   = s/1000
mss = [0.05, 0.5]
rss = [0.05, 0.1, 0.5, 1, 2, 4]
map(rss) do rs
    P = plot(title="\$r/s = $(round(rs, digits=2))\$")
    map(enumerate(mss)) do (i,ms)
        m = ms*s
        r = rs*s
        d = Fwd.distance(r)
        C = d*5
        xs = [C/2-d/2, C/2+d/2]
        q1, q2, D = TwoLocusModels.BA11(m, s, s, r)
        #AM = AeschbacherModel(m, [s for i=1:L], xs)
        #plot!(range(0, C, 500), x->1/Barriers.me(AM, x), 
        #    label="\$m/s=$ms, q=$(round(q1, digits=3)), D=$(round(D, digits=2))\$",
        #    yscale=:log10, color=i, lw=2, alpha=0.2)
        plot!(range(0, C, 500), x->scpred(m, xs, D, q1, q2, x), label="",
            yscale=:log10, color=i, ls=:dash, lw=1, 
            legend=:topright)
        # Predict allele frequencies
        loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M)
        Eq = 1 - EM.Ep[1]
        # ---------------------------
        #plot!(range(0, C, 500), 
        #    x->scpred(m, xs, 0, Eq, Eq, x), label="",
        #    yscale=:log10, color=i, lw=1,
        #    legend=:topright)
        #AM = AeschbacherModel(m, [s*EM.Ep[i] for i=1:L], xs)
        #plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="",
        #    yscale=:log10, color=i, lw=1, ls=:dot)
        BP = Equilibrium(BPModel(m=m, s=fill(s,L), u=u, Ne=float(NB), xs=xs))
        plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="",
            yscale=:log10, color=i, lw=1, ls=:dot)
        y1,y2 = ylims(P)
        plot!(ylim=(-Inf,y2*5))
    end 
    P
end |> x->plot(x..., size=(900,500), xlabel="map position (M)", ylabel="\$T\$", margin=3Plots.mm)

mss = [0.05, 0.5]
rss = [0.1, 0.5, 1, 2]
res = map(rss) do rs
    map(enumerate(mss)) do (i,ms)
        m   = ms*s
        r   = rs*s
        d   = Fwd.distance(r)
        C   = 5d
        xs  = [C/2-d/2, C/2+d/2]
        L   = 2
        R   = LinearMap(C)
        AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
        AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
        MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
        MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
        NA  = 1
        nA = collect(1:NA)
        nB = collect(1:NB) .+ NA
        xA = [ ones(Bool, L) for _=1:NA]
        xB = [zeros(Bool, L) for _=1:NB]
        popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
        popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
        AM = AeschbacherModel(m, [s for i=1:L], xs)
        TAM = 1/Barriers.me(AM, xs[1]*0.99)
        ngen = ceil(Int64, 2TAM)
        @info ngen
        res = pmap(1:10) do _
            mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
            mpop, ts, qs = simulate!(
                mpop, init_ts(mpop), ngen, x->twolocus_cb(x.popB, 1, 2), 
                every=ngen÷1000)
        end
    end
end

tres = map(zip(res, rss)) do (res_, rs)
    map(enumerate(zip(res_, mss))) do (i,(res__, ms))
        tbs = map(res__) do (_,ts,_)
            _ts = Fwd._add_grand_ancestor(ts)
            Fwd.diffdiv(_ts)[[1,4]]
        end
        xx, yy = Fwd.summarize_wins(tbs)
        yb = vec(mean(yy, dims=1))
        xx, yb
    end
end

#serialize("data/twolocus-2025-12-19.jls", (mss, rss, tres))
mss, rss, tres = deserialize("data/twolocus-2025-12-19.jls")

Ps = map(zip(tres, rss)) do (res_, rs)
    P = plot(title="\$r/s = $(round(rs, digits=2))\$")
    map(enumerate(zip(res_, mss))) do (i,(res__, ms))
        m   = ms*s
        r   = rs*s
        d   = Fwd.distance(r)
        C   = 5d
        xs  = [C/2-d/2, C/2+d/2]
        L   = 2
        # --------------------------------
        xx, yb = res__
        plot!(P, xx, yb, color=:black, label="")
        # estimate q by MC ---------------
        mod = TwoLocusModels.HaploidMainlandIsland(
            m=m, w=[1, 1-s, 1-s, (1-s)^2], c=r, u=u)
        rng = Random.default_rng()
        x = [NB, 0, 0, 0]
        xm = x
        nn = 1000000
        for _=1:nn
            x = generation(rng, mod, x)
            xm += x
        end
        q1, q2, D = TwoLocusModels.qd((xm ./ NB) ./ nn)
        @info q1, q2, D
        # --------------------------------
        #q1, q2, D = TwoLocusModels.BA11(m, s, s, r)
        #@info q1, q2, D
        # --------------------------------
        loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
        R = Fwd.rec_matrix(xs)
        A = Barriers.Architecture(loci, xs, R)
        M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
        EM = Barriers.Equilibrium(M)
        Eq = 1 - EM.Ep[1]
        # --------------------------------
        AM = AeschbacherModel(m, [s for i=1:L], xs)
        plot!(range(0, C, 500), x->1/Barriers.me(AM, x), label="",
            yscale=:log10, color=i, ls=:dot)
        ## --------------------------------
        BP = Equilibrium(BPModel(m=m, s=[s, s], xs=xs, Ne=NB, u=u))
        BP.Ep .= 1 .- [q1, q2]
        plot!(range(0, C, 500), x->1/Barriers.me(BP, x), label="",
            yscale=:log10, color=i, lw=1, ls=:dash)
        # --------------------------------
        plot!(range(0, C, 500), x->scpred(m, xs, D, q1, q2, x), label="",
            yscale=:log10, color=i, lw=5, alpha=0.3, 
            legend=:topright)
        # --------------------------------
        #plot!(range(0, C, 500), 
        #    x->scpred(m, xs, 0, Eq, Eq, x), 
        #    yscale=:log10, color=i, lw=1,
        #    label="\$m/s=$ms, q=$(round(q1, digits=3)), D=$(round(D, digits=2))\$",
        #    legend=:topright)
        # ---------------------------
        #plot!(range(0, C, 500), 
        #    x->scqle(m, xs, Eq, Eq, s, s, x), 
        #    yscale=:log10, color=i, lw=1,
        #    label="\$m/s=$ms, q=$(round(q1, digits=3)), D=$(round(D, digits=2))\$",
        #    legend=:topright)
        # ---------------------------
        y1,y2 = ylims(P)
        plot!(ylim=(-Inf,y2))
    end
    P
end 

annotate!(Ps[1], 0., 10^4/3, text("\$m/s=0.5\$", 8, :left))
annotate!(Ps[1], 0., 10^5/2, text("\$m/s=0.05\$", 8, :left))
plot(Ps..., size=(620,380), 
    xlabel="map position (M)", ylabel="\$T\$", margin=3Plots.mm)



res = map(1:10) do _
    NB  = 500
    s   = 10/NB
    u   = s/1000
    m   = s*0.2
    r   = s/2
    d   = Fwd.distance(r)
    C   = 5d
    xs  = [C/2-d/2, C/2+d/2]
    L   = 2
    R   = LinearMap(C)
    AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
    NA  = 1
    nA = collect(1:NA)
    nB = collect(1:NB) .+ NA
    xA = [ ones(Bool, L) for _=1:NA]
    xB = [zeros(Bool, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
    ngen = 10^6
    mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
    mpop, ts, qs = simulate!(
        mpop, init_ts(mpop), ngen, x->twolocus_cb(x.popB, 1, 2))
end

tss = map(getindex.(res, 2)) do ts
    TS.diffdiv(TS._add_grand_ancestor(ts))
end
xx, tbs = TS.summarize_wins(getindex.(tss, Ref([1,3])))
tb = mean(tbs, dims=1) |> vec

BP = Equilibrium(BPModel(m=m, s=[s,s], u=u, Ne=float(NB), xs=xs))
_tab(m) = 1/m + NA
_tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)

plot(xx, tb, yscale=:log10)
pred = map(range(0, C, 500)) do x
    r_1 = Fwd.recrate(abs(xs[1] - x))
    r_2 = Fwd.recrate(abs(xs[2] - x))
    p_1 = BP.Ep[1]; p_2 = BP.Ep[2]
    s_1 = s_2 = s
 #   xx = m*(p_1*s_1 + r_2)/(p_1^2*s_1^2 + p_1*p_2*s_1*s_2 + p_1*r_1*s_1 + p_1*r_2*s_1 + p_2*r_1*s_2 + r_1*r_2)
    me = Barriers.me(BP, x)
    xx = if x < xs[1]
        me/r_1
    elseif x > xs[2]
        me/r_2    
    else
        #me*(r_1 + r_2)/(r_1*r_2)
        me/(s_1*p_1 + s_2*p_2 + r_1 + r_2) +
            me*(p_1*s_1 + r_1 + r_2)*(p_2*s_2 + r_1 + r_2)/(r_1*r_2*(s_1*p_1 + s_2*p_2 + r_1 + r_2))
    end
    yy = me/min(r_1, r_2)
    zz = me/(1/(1/r_1 + 1/r_2)) 
    #@info _xx, xx
    tw1 = (1-2xx)*_tw(NA,NB,me) + 2xx*_tab(me)
    tw2 = (1-2xx)*_tw(NA,NB*(1-xx),me) + 2xx*_tab(me)
    tw3 = _tw(NA,NB*(1-xx),me)
    tw4 = _tw(NA,NB,me)
    tw5 = (1-2yy)*_tw(NA,NB,me) + 2yy*_tab(me)
    tw6 = (1-2zz)*_tw(NA,NB,me) + 2zz*_tab(me)
    x, xx, tw1, tw2, tw3, tw4, tw5, tw6
end
plot!(getindex.(pred, Ref([1,3])), lw=2)
plot!(getindex.(pred, Ref([1,5])), lw=2)
plot!(getindex.(pred, Ref([1,7])), lw=2)
plot!(getindex.(pred, Ref([1,8])), lw=2)
#plot!(getindex.(pred, Ref([1,3])), lw=2)

plot(getindex.(pred, Ref([1,2])), lw=2)


    
mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
mpop, qs = simulate!(mpop, 500_000, x->twolocus_cb(x.popB, 1, 2))

Q = permutedims(hcat(qs...))
Q1 = Q[:,3] .+ Q[:,4]
Q2 = Q[:,2] .+ Q[:,4]
DD = Q[:,1] .* Q[:,4] .- Q[:,2] .* Q[:,3]
Q = [Q hcat(Q1, Q2, DD)]

tbs = map(res) do (_,ts,_)
    _ts = TS._add_grand_ancestor(ts)
    # should make sure everything has coalesced...
    TS.diffdiv(_ts)[[1,4]]
end
xx, yy = TS.summarize_wins(tbs)
yb = vec(mean(yy, dims=1))
xx, yy, yb
#serialize("data/twolocus-2025-12-10.jls", (xx,yy))

xx, yy = deserialize("data/twolocus-2025-12-03.jls")
yb = vec(mean(yy, dims=1))
xx, yy, yb

Z = map(last.(res)) do zs
    Z = permutedims(hcat(zs...))
    D = Z[:,1] .* Z[:,4] .- Z[:,2] .* Z[:,3]
    q1 = Z[:,3] .+ Z[:,4]
    q2 = Z[:,2] .+ Z[:,4]
    hcat(Z, D, q1, q2)
end
mn = vec(mean(mean(Z), dims=1))

q1, q2, D = TwoLocusModels.BA11(m, s, s, r)

xspan = (0.0,C)
plot(xx, yb, xlim=xspan, linetype=:steppost, color=:gray, alpha=0.9, 
    xlabel="map position (M)", ylabel="\$T\$", label="", margin=3Plots.mm)
vline!([xs], legend=:topright, label="", size=(500,220))
hline!([ngen], color=:black, ls=:dot, label="")
#AM = AeschbacherModel(m, [s for i=1:L], xs)
#plot!(range(0, C, 500), x->1/Barriers.me(AM, x), 
#    yscale=:log10, label="AB14", ls=:solid, color=1, lw=4, alpha=0.2)
plot!(range(0, C, 500), x->scpred(m, xs, D, q1, q2, x), 
    yscale=:log10, label="SC det.", color=:black, lw=5, alpha=0.2)
plot!(range(0, C, 500), x->scpred(m, xs, mn[5], mn[6], mn[7], x), 
    yscale=:log10, label="SC MC", color=:black, lw=1)
BP = Equilibrium(BPModel(m=m, s=[s for i=1:L], xs=xs, Ne=NB, u=u))
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), 
    yscale=:log10, label="BP pred.", color=1, lw=1, )
BP.Ep .= 1 .- mn[6:7]
plot!(range(0, C, 500), x->1/Barriers.me(BP, x), 
    yscale=:log10, label="BP MC", color=1, lw=1, ls=:dash )
plot!(title="\$m/s = $(m/s), r/s=$(r/s), Ns=$(NB*s)\$")

loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
EM = Barriers.Equilibrium(M)
Eq = 1 - EM.Ep[1]
plot(xx, yb, xlim=xspan, linetype=:steppost, color=:gray, alpha=0.9, 
    xlabel="map position (M)", ylabel="\$T\$", label="", margin=3Plots.mm)
plot!(range(0, C, 500), x->scpred(m, xs, 0, Eq, Eq, x), 
    yscale=:log10, label="SC, \$\\mathbb{E}[q],D=0\$", color=3, lw=4, alpha=0.3)
plot!(range(0, C, 500), x->scpred(m, xs, D, q1, q2, x), 
    yscale=:log10, label="SC, \$\\tilde{q},\\tilde{D}\$", color=4, lw=4, alpha=0.3)
AM2 = AeschbacherModel(m, [s*EM.Ep[1] for i=1:L], xs)
plot!(range(0, C, 500), x->1/Barriers.me(AM2, x), 
    yscale=:log10, label="AB14", ls=:dot, color=:black, lw=1, )
vline!([xs], legend=:topright, label="", size=(500,220))
hline!([ngen], color=:black, ls=:dot, label="")
plot!(title="\$m/s = $(round(m/s, digits=2)), r/s=$(r/s), Ns=$(NB*s)\$")

function AB14full(r1, r2, s1, s2, q1, q2, D)
    A = s1*s2*D + (s1*q1 - r1)*(s2*q2 - r2)
    B = 2s1*s2*D + (2s1*q1 - s1 - r1)*(2s2*q2 - s2 - r2)
    A/B
end

plot(xx, yb, xlim=xspan, linetype=:steppost, color=:gray, alpha=0.9, 
    xlabel="map position (M)", ylabel="\$T\$", label="", margin=3Plots.mm)
AM2 = AeschbacherModel(m, [s for i=1:L], xs)
plot!(range(0, C, 500), x->1/Barriers.me(AM2, x), 
    yscale=:log10, label="AB14", color=:black, lw=1, )
plot!(range(xs[1]+1e-2, xs[2]-1e-2, 100), x->1/(m*AB14full(
    Fwd.recrate(x-xs[1]), Fwd.recrate(xs[2]-x), s, s, q1, q2, D)),
    yscale=:log10, label="AB14", ls=:dot, color=:black, lw=1, )


# This is for the case with the neutral locus on the left
function twolocus_rv(s1, s2, r, r1, q1, q2, D, kmax=100) 
    p1 = 1-q1
    p2 = 1-q2
    p1k = 0.0 
    p2k = 0.0
    Dk = 0.0
    wres = exp(-s1*q1 - s2*q2)
    wk(p1, p2, D) = (p1^2 + D) + (p1*(1-p1) - D)*exp(-s2) + ((1-p1)*p2 - D)*exp(-s1) + ((1-p1)*(1-p2) + D)*exp(-s1 - s2) 
    g = wk(p1k, p2k, Dk)/wres
    for k=1:kmax
        Dk = (1-r)*Dk + r1*(D - Dk) + r1*(1-r-r1)*(p1 - p1k)*(p2 - p2k)
        p1k = p1k + r1*(p1 - p1k)
        p2k = p2k + (r+r1)*(p2 - p2k)
        g *= wk(p1k, p2k, Dk)/wres
    end
    return g
end

function twolocus_rvx(x, xs, s1, s2, q1, q2, D)
    if x < xs[1] || x > xs[2]
        r1 = Fwd.recrate(min(abs(xs[1] - x), abs(xs[2]-x)))
        r = Fwd.recrate(xs[2]-xs[1])
        twolocus_rv(s1, s2, r, r1, q1, q2, D)
    else
        NaN
    end
end

twolocus_rv(s, s, r, r, q1, q2, D)

loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
R = Fwd.rec_matrix(xs)
A = Barriers.Architecture(loci, xs, R)
M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
EM = Barriers.Equilibrium(M)
Eq = 1 - EM.Ep[1]
plot(xx, yb, xlim=xspan, linetype=:steppost, color=:gray, alpha=0.9, 
    xlabel="map position (M)", ylabel="\$T\$", label="", margin=3Plots.mm)
plot!(range(0, C, 500), x->scpred(m, xs, D, Eq, Eq, x), 
    yscale=:log10, label="SC, \$\\mathbb{E}[q],D=\\hat{D}\$", color=4, lw=4, alpha=0.3)
plot!(range(0, C, 500), x->1/(m*twolocus_rvx(x, xs, s, s, q1, q2, D)), 
    yscale=:log10, label="SC, \$\\mathbb{E}[q],D=\\hat{D}\$", color=1, lw=2, alpha=0.3)
plot!(range(0, C, 500), x->1/Barriers.me(EM, x), 
    yscale=:log10, label="SC, \$\\mathbb{E}[q],D=\\hat{D}\$", color=:black, lw=2, alpha=0.3)

# ---------------------------------------------------------------
# Vary mig. rates
#ms = [s/50, s/20, s/10, s/5]
#ms = [s/5]
#nrep = 1
#res = map(ms) do m
#    res = map(1:nrep) do _
#        mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
#        mpop, ts, qs = simulate!(
#            rng, mpop, init_ts(mpop), ngen, x->mean(x.popB.x))
#    end
#end
#
#q1 = mean(res[1][1][end])
#
#T = map(res) do X
#    tbs = map(X) do (_,ts,_)
#        _ts = Fwd._add_grand_ancestor(ts)
#        # should make sure everything has coalesced...
#        Fwd.diffdiv(_ts)[[1,4]]
#    end
#    xx, yy = Fwd.summarize_wins(tbs)
#    yb = vec(mean(yy, dims=1))
#    xx, yy, yb
#end
#
#plot()
#map(zip(ms, T)) do (m, t)
#    plot!(t[1], t[3], linetype=:steppost, label="\$m=$m\$")
#end
#vline!([xs], legend=:outertopright, label="", size=(700,200))
#hline!([ngen])
#AM = AeschbacherModel(ms[1], [s for i=1:L], xs)
#plot!(range(0, C, 200), x->1/Barriers.me(AM, x), yscale=:log10)
#
#plot(xx, 1 ./ yb)
#
#T = map(res) do X
#    tbs = map(X) do (_,ts,_)
#        Fwd.theights(ts)
#    end
#    xx, yy = Fwd.summarize_wins(tbs)
#    yb = vec(mean(yy, dims=1))
#    xx, yy, yb
#end
#
#plot()
#map(enumerate(zip(ms, T))) do (i,(m, t))
#    plot!(t[1], t[3], linetype=:steppost, label="\$m=$m\$", color=i)
#    AM = AeschbacherModel(m, [-s for i=1:L], xs)
#    plot!(range(0, C, 200), x->1/Barriers.me(AM, x), 
#        label=false,
#        color=i, ls=:dash, yscale=:log10)
#end
#vline!([xs], legend=:outertopright, label="", size=(700,200))
#hline!([ngen])
#
#
## self-consistent solution?
#gff12(r, q, m) = r*q/((1-q)*m + r*q)
#
#function predict_twolocus(m, r, s, N, u)
#    d(m) = Wright(2N*s, N*(m+u), N*u, 0.5)
#    q = mean(d(m))
#    while true
#        g = gff12(r, q, m)
#        _q = mean(d(m*g))
#        abs(q - _q) < 1e-5 && return _q
#        q = _q
#    end
#end
#
## check how this approach fits for increasing r/s
## note that it cannot work in general: it predicts zero mₑ when fully linked
## to a barrier, and then using mₑ/s to predict p must break down
#
#ngen = 500NB
#rng = Random.seed!(12)
#rss = exp2.(range(-2.5, 2.5, 10))
#m = 0.2*s
#res = map(rss) do rs
#    r  = min(0.5-1e-2,s*rs)
#    d  = Fwd.distance(r)
#    @info r, d
#    xs = [0., d]
#    R  = LinearMap(d)
#    AA = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
#    AB = Architecture([BiAllelic(u) for _=1:L], xs, R)
#    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
#    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
#    NA  = 1
#    nA = collect(1:NA)
#    nB = collect(1:NB) .+ NA
#    xA = [ ones(Bool, L) for _=1:NA]
#    xB = [zeros(Bool, L) for _=1:NB]
#    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
#    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
#    mpop = Fwd.TwoPopOneWay(m, popA, popB)
#    mpop, ts, qs = simulate!(
#        rng, mpop, init_ts(mpop), ngen, x->mean(x.popB.x))
#end
#
#rss2 = exp2.(range(-2.5, 2.5, 200))
#res2 = map(rss2) do rs
#    r  = min(0.5-1e-2,s*rs)
#    qp = predict_twolocus(m, r, s, NB, u)
#end
#
#ps = map(x->mean(mean(x)), last.(res)) 
#plot(rss2, res2, xscale=:log2, color=:black, ylabel="\$q\$", label="iterative SC")
#scatter!(rss, ps, xscale=:log2, xlabel="\$r/s\$", 
#    label="simulation", legend=:bottomleft,
#    title="\$N = $NB, s=$s, m=$(round(m,digits=2)), u=$u\$")
#
#d = Wright(2NB*s, NB*(m+u), NB*u, 0.5)
#plot(0:0.01:1, x->pdf(d, x))
#vline!([mean(d)])
#vline!([m/s])
#
#_, ts, qs = res[end]
