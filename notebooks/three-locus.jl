
using Distributed; addprocs(10)
@everywhere using Fwd, TreeSequences, StatsBase
@everywhere using WrightDistribution, Barriers
using Plots; plotsdefault()
using Serialization
import TreeSequences as TS

@everywhere const states=[[0,0,0],[0,0,1],[0,1,0],[1,0,0],[0,1,1],[1,0,1],[1,1,0],[1,1,1]]

@everywhere function threelocus_cb(pop)
    pm = proportionmap(pop.x)
    [haskey(pm, x) ? pm[x] : 0.0 for x in states]
end
    
@everywhere function threelocus_model(N, r, m, s, u=s/500)
    d   = Fwd.distance(r)
    L   = 3
    flank = L*d
    xs  = flank .+ [k*d for k=0:L-1]
    C   = last(xs) + flank
    R   = LinearMap(C)
    AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
    NA  = 1
    nA = collect(1:NA)
    nB = collect(1:N) .+ NA
    xA = [ ones(Bool, L) for _=1:NA]
    xB = [zeros(Bool, L) for _=1:N]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=N, arch=AB, gpm=MB, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
end

@everywhere function threelocus_stats(x)
    z = zeros(3)
    for k=1:3
        z[k] = sum(x[filter(i->states[i][k] == 0, 1:length(states))])
    end
    # add LD...
    z 
end

function _threelocus_prediction(p, N, m, r, s, u)
    p_1, p_2, p_3 = p
    r_1, r_2 = r
    s_1, s_2, s_3 = s
    g1 = (p_1*s_1 + r_1)*(p_1*s_1 + p_2*s_2 + r_2 + r_1)/((p_1*s_1 + p_2*s_2 + r_1)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_2 + r_1))
    g2 = (p_1*p_2*p_3*s_1*s_2*s_3 + (p_2*s_2 + r_1)*(p_2*s_2 + r_2)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))/((p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
    g3 = (p_3*s_3 + r_2)*(p_3*s_3 + p_2*s_2 + r_2 + r_1)/((p_3*s_3 + p_2*s_2 + r_1)*(p_3*s_3 + p_2*s_2 + p_1*s_1 + r_2 + r_1))
    d1 = Wright(2N*s_1, N*(u + g1*m), N*u, 0.5)
    d2 = Wright(2N*s_2, N*(u + g2*m), N*u, 0.5)
    d3 = Wright(2N*s_3, N*(u + g3*m), N*u, 0.5)
    p = 1 .- [mean(d1), mean(d2), mean(d3)]
end

function threelocus_prediction(N, m, r, s, u, tol=1e-5, α=0.2)
    p = ones(3)
    while true 
        p_ = _threelocus_prediction(p, N, m, r, s, u)
        sum(abs.(p .- p_)) < tol && return p_
        p = p*(1-α) .+ p_ * α
    end
end

let s=0.02, N=10/s, rss=[0.1, 0.2, 0.5, 1.0, 5.0], u=s/200, mss=range(0.01, 1.8, 100)
    Ps = [plot(), plot()]
    map(enumerate(rss)) do (i,rs)
        d = Fwd.distance(rs*s)
        r = [rs*s, rs*s]
        xs = [0, d, 2d]
        ps = map(mss) do ms
            BP = BPModel(Ne=N, s=[s,s,s], u=u, m=ms*s, xs=xs)
            Eq = Equilibrium(BP, α=0.2).Ep
        end 
        plot!(Ps[1], mss, pcat(ps...)[:,1], color=i, ls=:dash, label="")
        plot!(Ps[2], mss, pcat(ps...)[:,2], color=i, ls=:dash, label="")
        ps = map(mss) do ms
            BP = BPModel(Ne=N, s=[s,s,s], u=u, m=ms*s, xs=xs)
            Eq = Barriers.EquilibriumS(BP, α=0.2).Ep
        end 
        plot!(Ps[1], mss, pcat(ps...)[:,1], color=i, ls=:dot, label="")
        plot!(Ps[2], mss, pcat(ps...)[:,2], color=i, ls=:dot, label="")
        ps = map(ms->threelocus_prediction(N, ms*s, r, [s,s,s], u), mss)
        plot!(Ps[1], mss, pcat(ps...)[:,1], color=i, label="\$r/s=$rs\$")
        plot!(Ps[2], mss, pcat(ps...)[:,2], color=i, label="\$r/s=$rs\$")
    end
    plot!(Ps..., legend=:topright, size=(550,230), 
        xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
end

res = let s=0.02, Ns=10, rss=[0.1, 0.2, 0.5, 1.0, 5.0], 
    u=s/200, mss=range(0.01, 1.8, 15)
    map(rss) do rs
        @info rs
        qs = pmap(mss) do ms
            m = ms*s
            N = ceil(Int, Ns/s)
            model = threelocus_model(N, rs*s, m, s, u)
            ngen = 5000N
            evry = ceil(Int,ngen÷2500)
            mpop, qs = simulate!(model, ngen, x->threelocus_cb(x.popB), every=evry)
            pcat(map(threelocus_stats, qs)...)
        end
        (Ns=Ns, rs=rs, ms=mss, qs=qs)
    end
end 

let s=0.02, N=10/s, rss=[0.1, 0.2, 0.5, 1.0, 5.0], u=s/200, mss=range(0.01, 1.8, 100)
    Ps = [plot(), plot()]
    map(enumerate(rss)) do (i,rs)
        d = Fwd.distance(rs*s)
        r = [rs*s, rs*s]
        xs = [0, d, 2d]
        ps = map(mss) do ms
            BP = BPModel(Ne=N, s=[s,s,s], u=u, m=ms*s, xs=xs)
            Eq = Equilibrium(BP, α=0.2).Ep
        end 
        plot!(Ps[1], mss, pcat(ps...)[:,1], color=i, ls=:solid, label="\$r/s=$rs\$")
        plot!(Ps[2], mss, pcat(ps...)[:,2], color=i, ls=:solid, label="\$r/s=$rs\$")
        ps = map(ms->threelocus_prediction(N, ms*s, r, [s,s,s], u), mss)
        plot!(Ps[1], mss, pcat(ps...)[:,1], color=i, ls=:dash, label="")
        plot!(Ps[2], mss, pcat(ps...)[:,2], color=i, ls=:dash, label="")
        se = [mcse(Chains(res[i].qs[k])).nt.mcse[1:2] for k=1:length(res[i].qs)]
        scatter!(Ps[1], res[i].ms, map(x->mean(x[:,1]), res[i].qs), err = 2first.(se), color=i, ms=3, msc=i, lw=1, label="")
        scatter!(Ps[2], res[i].ms, map(x->mean(x[:,2]), res[i].qs), err = 2last.(se), color=i, ms=3, msc=i, lw=1, label="")
    end
    plot!(Ps..., legend=:topright, size=(550,230), 
        xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
end

res = let s=0.02, Ns=10, rs=0.1, u=s/200, ms=1.3
    @info rs
    m = ms*s
    N = ceil(Int, Ns/s)
    model = threelocus_model(N, rs*s, m, s, u)
    ngen = 5000N
    evry = ceil(Int,ngen÷2500)
    mpop, qs = simulate!(model, ngen, x->threelocus_cb(x.popB), every=evry)
    pcat(map(threelocus_stats, qs)...)
    (Ns=Ns, rs=rs, ms=mss, qs=qs)
end 

Eq = let s=0.02, Ns=10, rs=0.1, u=s/200, ms=1.3
    d = Fwd.distance(rs*s)
    r = [rs*s, rs*s]
    xs = [0, d, 2d]
    BP = BPModel(Ne=Ns/s, s=[s,s,s], u=u, m=ms*s, xs=xs)
    Eq = Equilibrium(BP, α=0.2).Ep
end


# ==================================================================================
mpop, ts, qs = simulate!(
    mpop, init_ts(mpop), ngen, x->threelocus_cb(x.popB), every=ngen÷5000)

plot(TS.diffdiv(TS._add_grand_ancestor(ts))[[1,4]], yscale=:log10)


P = pcat(map(threelocus_stats, qs)...)
plot(P)

# BP estimate for middle locus del. allele frequency
function q3(p, s, r_1, r_2, m)
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    m * (p_1*p_2^2*s_1*s_2^2 + p_1*p_2*p_3*s_1*s_2*s_3 + p_1*p_2*r_1*s_1*s_2 + p_1*p_2*r_2*s_1*s_2 + p_1*r_1*r_2*s_1 + p_2^3*s_2^3 + p_2^2*p_3*s_2^2*s_3 + 2*p_2^2*r_1*s_2^2 + 2*p_2^2*r_2*s_2^2 + p_2*p_3*r_1*s_2*s_3 + p_2*p_3*r_2*s_2*s_3 + p_2*r_1^2*s_2 + 3*p_2*r_1*r_2*s_2 + p_2*r_2^2*s_2 + p_3*r_1*r_2*s_3 + r_1^2*r_2 + r_1*r_2^2)/(p_2*s_2*(p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
end

p = vec(mean(P,dims=1))
q3(p, [s, s, s], r, r, m)


function g3(p, s, r_1, r_2) 
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_1*p_2^2*s_1*s_2^2 + p_1*p_2*p_3*s_1*s_2*s_3 + p_1*p_2*r_1*s_1*s_2 + p_1*p_2*r_2*s_1*s_2 + p_1*r_1*r_2*s_1 + p_2^3*s_2^3 + p_2^2*p_3*s_2^2*s_3 + 2*p_2^2*r_1*s_2^2 + 2*p_2^2*r_2*s_2^2 + p_2*p_3*r_1*s_2*s_3 + p_2*p_3*r_2*s_2*s_3 + p_2*r_1^2*s_2 + 3*p_2*r_1*r_2*s_2 + p_2*r_2^2*s_2 + p_3*r_1*r_2*s_3 + r_1^2*r_2 + r_1*r_2^2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
end

function g3b(p, s, r_1, r_2)
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_2*s_2 + r_1)*(p_2*s_2 + r_2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2))
end

function g4(p, s, r_1, r_2)
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_1*s_1 + r_1)*(p_1*s_1 + p_2*s_2 + r_1 + r_2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
end

gs1 = map(p-> g3(p, [s,s,s], r, r), eachrow(P))
gs2 = map(p-> g3b(p, [s,s,s], r, r), eachrow(P))
gs3 = map(p-> g4(p, [s,s,s], r, r), eachrow(P))

scatter(gs1, gs2)
plot!(x->x)

d1 = Wright(2NB*s, NB*(u + mean(gs1)*m), NB*u, 0.5)
d2 = Wright(2NB*s, NB*(u + mean(gs2)*m), NB*u, 0.5)
d3 = Wright(2NB*s, NB*(u + mean(gs3)*m), NB*u, 0.5)
stephist(1 .- P, bins=50, norm=true, label="", legend=:topright)
plot!(range(0,1,200), x->pdf(d1, x), label="derived")
plot!(range(0,1,200), x->pdf(d2, x), label="guess")
plot!(range(0,1,200), x->pdf(d3, x), label="left")


# probably Nm more important than Ns, as long as Ns appreciable
res = let s=0.02, 
    Nss=[5,10,20], 
    rss=[0.1, 1, 10], 
    mss=range(0, 1.5, 20)
    map(Nss) do Ns
        N = ceil(Int, Ns/s)
        ngen = 10N
        pmap(rss) do rs
            r = rs*s
            map(mss) do ms
                m = ms*s
                model = threelocus_model(N, r, m, s)
                mpop, ts, qs = simulate!(
                    model, init_ts(model), 
                    ngen, x->threelocus_cb(x.popB), 
                    every=ngen÷1000)
                (s=s, Ns=Ns, rs=rs, ms=ms, ts=ts, qs=qs)
            end
        end 
    end
end

map(res) do Y
    P = plot(title="\$Ns=$(Y[1][1].Ns)\$")
    map(enumerate(Y)) do (i,Xs)
        zs = map(Xs) do x
            @unpack Ns, ms, s, rs, qs = x
            p = mean(map(threelocus_stats, qs))
            s = 0.02
            q = q3(p, [s,s,s], rs*s, rs*s, ms*s)
            g = g3(p, [s,s,s], rs*s, rs*s)
            @info Ns * ms
            q > 1 && @warn "nonsense ($rs, $ms)" 
            1-p[2], q
        end 
        scatter!(P, zs, color=i, label="") #label="\$r/s=$(Xs[1].rs)\$", ms=3)
    end 
    plot!(P, x->x, color=:gray, legend=:topleft, label="")
    hline!([1], color=:black, ls=:dot, label="")
end |> x-> plot!(x..., size=(700,220), layout=(1,3), 
    xscale=:log10, yscale=:log10)

map(res) do Y
    P = plot(title="\$Ns=$(Y[1][1].Ns)\$")
    G = map(enumerate(Y)) do (i,Xs)
        gs = map(Xs) do x
            @unpack Ns, ms, s, rs, qs = x
            p = mean(map(threelocus_stats, qs))
            s = 0.02
            g = g3(p, [s,s,s], rs*s, rs*s)
            g
        end 
    end |> x->hcat(x...)
    rss = [y[1].rs for y in Y]
    plot!(P, rss, G', marker=true, ms=1.5, xscale=:log10, ylim=(0.3,1))
end |> x-> plot!(x..., size=(700,220), layout=(1,3), xlabel="\$r/s\$")

# the weird stuff is due to swamping.

function _getxs(r, L)
    d  = Fwd.distance(r)
    xs = [k*d for k=0:L-1]
end

map(res) do Y
    P = plot(title="\$Ns=$(Y[1][1].Ns)\$")
    G = map(enumerate(Y)) do (i,Xs)
        gs = map(Xs) do x
            @unpack Ns, ms, s, rs, qs = x
            p = mean(map(threelocus_stats, qs))
            s = 0.02
            g = g3(p, [s,s,s], rs*s, rs*s)
            q = q3(p, [s,s,s], rs*s, rs*s, ms*s)
            pp = Equilibrium(BPModel(
                m=ms*s, s=[s,s,s], xs=_getxs(rs*s,3), 
                Ne=Ns/s, u=s/500), nmax=101).Ep[2]
            ms, p[2], max(1-q, 0), pp
            # add mₑ based self-consistent prediction, compare with
            # 'empirical' mₑ based prediction
        end
        plot!(first.(gs), getindex.(gs,2), 
            color=i, label="\$r/s=$(Xs[1].rs)\$", legend=:topright)
        scatter!(first.(gs), getindex.(gs,3), 
            color=i, label="", legend=:topright, ms=2)
        scatter!(first.(gs), getindex.(gs,4), 
            color=i, label="", marker=:x, legend=:topright, ms=2)
    end 
    P
end |> x-> plot!(x..., size=(700,220), layout=(1,3), xlabel="\$m/s\$", 
    margin=3Plots.mm)


# -----------------
s = 0.02
u = s/200
Nss = [5, 10]
rss = [0.1, 1., 10.]
mss = Dict(Nss[1]=>(0.01, 1.2), Nss[2]=>(0.01,1.5))

res = map(Nss) do Ns
    map(rss) do rs
        qs = pmap(range(mss[Ns]..., 15)) do ms
            m = ms*s
            N = ceil(Int, Ns/s)
            model = threelocus_model(N, rs*s, m, s, u)
            ngen = 5000N
            evry = ceil(Int,N÷5)
            mpop, qs = simulate!(model, ngen, x->threelocus_cb(x.popB), every=evry)
            ms, qs
        end
        (Ns=Ns, rs=rs, qs=qs)
    end
end 

res2 = map(Nss) do Ns
    map(rss) do rs
        ps = map(range(mss[Ns]..., 100)) do ms
            m = ms*s
            p = Equilibrium(BPModel(
                m=ms*s, s=[s,s,s], xs=[0, rs*s, 2rs*s], 
                Ne=Ns/s, u=u), α=0.2, nmax=1000).Ep[2]
            ms, p
        end
        (Ns=Ns, rs=rs, ps=ps)
    end
end 

PS = map(1:length(res)) do k
    P = plot()
    G = map(1:length(res[k])) do i
        plot!(res2[k][i][3], color=i)
        xx = first.(res[k][i][3])
        yy = map(threelocus_stats ∘ mean, last.(res[k][i][3]))
        scatter!(xx, getindex.(yy, 2), color=i, ms=2)
    end
    P
end 
plot(PS..., layout=(1,2), size=(500,200))


plot(pcat(map(threelocus_stats, res[1][1].qs[12][2])...))

model = threelocus_model(500, 0.1*0.02, 1.1*0.02, 0.02)
mpop, qs = simulate!(model, 5_000_000, x->threelocus_cb(x.popB), every=200)

yy = mapreduce(threelocus_stats, hcat, qs)
plot(yy')

plot(PS[2])
scatter!([1.1], [mean(yy[2,:])], marker=:x, color=1)

            
# Numerical issues
Ns = 5.0
s  = 0.02 
Ne = Ns/s
m  = 0.7*s
u  = s/500
r  = 0.1*s
xs = [0, r, 2r]
BP = BPModel(m=m, s=[s,s,s], xs=xs, Ne=Ne, u=u)
Q  = Barriers._fixed_point_iteration(BP, ones(3), 1e-6, 30)
plot(Q)
# The oscillation is due to: left locus has low p -> middle gets higher g,
# but middle had high p, so left gets low g, next iteration vice versa...

BP = BPModel(m=m, s=[s,s,s], xs=xs, Ne=Ne, u=u)
Q  = Barriers._fixed_point_iteration(BP, ones(3), 1e-6, 1000, 0.1)
plot(Q)




