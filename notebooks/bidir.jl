using Distributed; addprocs(10)
@everywhere using Pkg; @everywhere Pkg.activate("/home/arzwa/dev/Fwd")
@everywhere using Fwd, Distributions, Random, Plots, Parameters, Serialization
import TreeSequences as TS
using Barriers
using PyCall
pg = pyimport("phasegen")

# Architecture
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
Ns   = 5. 
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
NB  = ceil(Int64, Ns/s̄)
NA  = NB
C   = 0.5
G   = C*100*10^6   # (1cM/Mb)
xs  = ceil.(Int64, ys .* G)
R   = LinearPhysMap(C=C, G=G)
    
# Set of nrep simulations
ms = 0.5
mAB = ms*s̄
mBA = ms*s̄
nrep = 10
res = pmap(1:nrep) do _
    arch = Architecture([BiAllelic(u) for _=1:L], xs, R)
    nA = collect(1:NA)
    nB  = collect(1:NB) .+ NA
    xA  = [ ones(Int, L) for _=1:NA]
    xB  = [zeros(Int, L) for _=1:NB]
    gpmA = GPMap([HaploidLocus( ss[i], i) for i=1:L])
    gpmB = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
    popA = WFPopulation(ploidy=Haploid(), 
        N=NA, arch=arch, gpm=gpmA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), 
        N=NB, arch=arch, gpm=gpmB, x=xB, nodes=nB)
    mpop = Fwd.MetaPop([popA, popB], [0. mAB; mBA 0.])
    mpop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), 100_000,
        pop->[mean(p.x) for p in pop.P], simplify=100, every=200)
end
    
#serialize("data/bidir/2026-03-10.1.jls", res)

BP = Equilibrium(BPModel(m=mAB, xs=ys*C, Ne=NB, s=ss, u=u))

BP2 = Barriers.BPTwoPop(
    m12 = mAB,
    m21 = mBA,
    s1  = ss,
    s2  = ss,
    xs  = ys * C,
    Ne1 = 1.0NA,
    Ne2 = 1.0NB,
    u   = u)

PP = Equilibrium(BP2, α=0.1, tol=1e-7)

Qs = map(res) do (_, _, qs)
    p1 = first.(qs) |> mean
    p2 = last.(qs) |> mean
    p1, p2
end
(q1, q2) = (mean(first.(Qs)), mean(last.(Qs)))
tss = getindex.(res, 2)

div = q1 .- q2

pa = PP.Ep[:,1]
pb = PP.Ep[:,2]
P1 = scatter(q1, pa, label="pop. A", xlabel="\$q\$ (sim.)", ylabel="\$q\$ (pred.)")
scatter!(q2, 1 .- pb, label="pop. B", xlim=(0,1), ylim=(0,1))
plot!(x->x, ls=:dot, color=:gray, label="")
title!(P1, @sprintf("\$Nm_{AB} = %.1f, Nm_{BA} = %.1f, L=%d, L\\bar{s}=%.2f, N=%d\$", mAB*NB, mBA*NA, L, Ls, NB))
P2 = scatter(pa, 1 .- pb, xlim=(0,1), ylim=(0,1), label="prediction")
scatter!(q1, q2, xlim=(0,1), ylim=(0,1), label="simulation", 
    xlabel="\$q_A\$", ylabel="\$q_B\$")
plot(P1, P2, legend=:topleft, size=(520,250), ms=2)

function bidir_fst(N1, N2, m21, m12)
    N = N1
    demography = pg.Demography(
        pop_sizes=Dict("A"=>N1/N, "B"=>N2/N),
        migration_rates=Dict(
            ("A","B")=>m21*N, 
            ("B","A")=>m12*N))
    coal = pg.Coalescent(
        n = Dict("A"=>1, "B"=>1),
        demography = demography)
    tab = coal.tree_height.mean * N
    coal = pg.Coalescent(
        n = Dict("A"=>2, "B"=>0),
        demography = demography)
    ta = coal.tree_height.mean * N
    coal = pg.Coalescent(
        n = Dict("B"=>2, "A"=>0),
        demography = demography)
    tb = coal.tree_height.mean * N
    Fst = 1 - ((ta + tb) / 2) / tab
    Fst, tab, ta, tb
end

# Fst estimates from simulatiosn with intervals
function estimate_fst(tss, ci=0.95)
    q0 = (1-ci)/2
    q1 = 1-q0
    tbs = map(tss) do ts
       _ts = TS._add_grand_ancestor(ts)
       TS.diffdiv(_ts)
    end
    ts = map(2:4) do k
        xx, yy = TS.summarize_wins(getindex.(tbs, [[1,k]]))
    end
    fst = 1 .- 0.5 .* (ts[1][2] .+ ts[2][2]) ./ ts[3][2]
    est = map(eachcol(fst)) do col
        mn = mean(col)
        q1 = quantile(col, (1-ci)/2)
        q2 = quantile(col, ci + (1-ci)/2)
        mn, mn - q1, q2 - mn
    end
    ts[1][1], getindex.(est, 1), getindex.(est, 2), getindex.(est, 3)
end

a, b, c, d = estimate_fst(tss, 0.9)

yy = map(range(0, C, 500)) do x
    y = G*x/C
    me21, me12 = Barriers.me(PP, x)
    y, bidir_fst(NA, NB, me21, me12)
end

plot(a, b, ribbon=(c,d), color=:gray, fillalpha=0.5, size=(700,200))
plot!(first.(yy), first.(last.(yy)), lw=1, color=:black)
sticks!(xs, div, lw=3, color=:firebrick, alpha=0.3)
plot!(xlabel="map position", ylabel="\$F_\\mathrm{ST}\$", margin=5Plots.mm)



# --------------------------
# m/s range


# Check theoretical predictions for some m/s range
mss = range(0.05, 5, 25)
preds = map(mss) do ms
    @info ms
    BP2 = Barriers.BPTwoPop(
        m12 = ms*s̄,
        m21 = ms*s̄,
        s1  = ss,
        s2  = ss,
        xs  = ys * C,
        Ne1 = 1.0NA,
        Ne2 = 1.0NB,
        u   = u)
    PP = Equilibrium(BP2, α=0.1, tol=1e-7)
    PP.Ep
end

map(1:2:L) do k
    plot(mss, mapreduce(p->p[k,:], hcat, preds)', ylim=(0,1))
    hline!([0.5])
end |> x->plot(x..., size=(800,800))

# Simulation
mss = range(0.5, 4, 10)
ress = pmap(enumerate(mss)) do (k, ms)
    @info ms
    nrep = 5
    res = map(1:nrep) do _
        arch = Architecture([BiAllelic(u) for _=1:L], xs, R)
        nA = collect(1:NA)
        nB  = collect(1:NB) .+ NA
        xA  = [ ones(Int, L) for _=1:NA]
        xB  = [zeros(Int, L) for _=1:NB]
        gpmA = GPMap([HaploidLocus( ss[i], i) for i=1:L])
        gpmB = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
        popA = WFPopulation(ploidy=Haploid(), 
            N=NA, arch=arch, gpm=gpmA, x=xA, nodes=nA)
        popB = WFPopulation(ploidy=Haploid(), 
            N=NB, arch=arch, gpm=gpmB, x=xB, nodes=nB)
        mpop = Fwd.MetaPop([popA, popB], [0. ms*s̄; ms*s̄ 0.])
        mpop, ts, qs = Fwd.simulate!(mpop, init_ts(mpop), 100_000,
            pop->[mean(p.x) for p in pop.P], simplify=100, every=200)
    end
    serialize("data/2026-03-10.$k.jls", res)
end
