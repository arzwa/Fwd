using Distributed#; addprocs(5)
@everywhere using Pkg; @everywhere Pkg.activate("/home/arzwa/dev/Fwd")
@everywhere using Fwd, Distributions, Random, Plots, Parameters, Serialization
import TreeSequences as TS
using Barriers

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
ms = 1.0
mAB = ms*s̄
mBA = ms*s̄
nrep = 5

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

res = deserialize("data/bidir/2026-03-10.1.jls")

mpop = res[1][1]
M = mpop.M

BP = Equilibrium(BPModel(m=mAB, xs=ys*C, Ne=NB, s=ss, u=u))

BP2 = Barriers.BPTwoPop(
    m12 = M[1,2],
    m21 = M[2,1],
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

scatter(div, pb .- (1 .- pa))
plot!(x->x, color=:lightgray, alpha=0.8, xlim=(0,1), ylim=(0,1))

function bidir_ct(N_A, N_B, m_AB, m_BA)
    ta = (2*N_A*N_B*m_AB^2 + 4*N_A*N_B*m_AB*m_BA + 2*N_A*N_B*m_BA^2 + N_A*m_AB + 3*N_A*m_BA)/(2*N_A*m_BA^2 + 2*N_B*m_AB^2 + m_AB + m_BA)
    tb = (2*N_A*N_B*m_AB^2 + 4*N_A*N_B*m_AB*m_BA + 2*N_A*N_B*m_BA^2 + 3*N_B*m_AB + N_B*m_BA)/(2*N_A*m_BA^2 + 2*N_B*m_AB^2 + m_AB + m_BA)
    tab = (2*N_A*N_B*m_AB^2 + 4*N_A*N_B*m_AB*m_BA + 2*N_A*N_B*m_BA^2 + N_A*m_AB + 2*N_A*m_BA + 2*N_B*m_AB + N_B*m_BA + 1)/(2*N_A*m_BA^2 + 2*N_B*m_AB^2 + m_AB + m_BA)
    fst = 1 - 0.5*(ta + tb)/tab
    fst, ta, tb, tab
end

# Fst estimates from simulations with intervals
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
    ts, ts[1][1], getindex.(est, 1), getindex.(est, 2), getindex.(est, 3)
end


function rho(BP, kmax=floor(Int, log2(0.5/BP.model.m12)))
    @unpack m12, m21, s1, s2, xs = BP.model
    Δ = BP.Ep[:,2] .- (1 .- BP.Ep[:,1])
    ws1 = map(0:kmax-1) do k
        exp(-sum([s1[i]*Δ[i] for i=1:length(xs)])/2^k)
    end 
    ws2 = map(0:kmax-1) do k
        exp(-sum([s2[i]*Δ[i] for i=1:length(xs)])/2^k)
    end 
    fs1 = m21*cumprod(2ws1)
    fs2 = m12*cumprod(2ws2)
    ρ1 = sum([1/2^k * fs1[k] for k=1:kmax])
    ρ2 = sum([1/2^k * fs2[k] for k=1:kmax])
    ρ1, ρ2
end


# Coalescence times/Fst
Ts = map([2,3,4]) do k
    a, b, c, d = estimate_coaltimes(tss, idx=k)
end
push!(Ts, estimate_fst(tss, 0.90)[2:end])

ρ1, ρ2 = rho(PP)

# between pop T
yy = map(range(0, C, 500)) do x
    y = G*x/C
    me21, me12 = Barriers.me(PP, x)
    _, ta, tb, tab = bidir_ct(NA, NB, me21, me12)
    y, ρ1*tb + ρ2*ta + (1-ρ1-ρ2)*tab, tab
end
plot(Ts[3][1:2], ribbon=(Ts[3][3:4]...,), size=(800,200), color=:lightgray)
plot!(getindex.(yy, Ref([1,2])), yscale=:log10, lw=2)
plot!(getindex.(yy, Ref([1,3])), yscale=:log10, lw=2)
# difference between predictions is minute

# within pop T
tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
yy = map(range(0, C, 500)) do x
    y = G*x/C
    me21, me12 = Barriers.me(PP, x)
    _, ta, tb, tab = bidir_ct(NA, NB, me21, me12)
    tbw = tw(NA, NB, me12)
    tb1 = 2ρ2*tab + (1-2ρ2)*tb
    tb2 = 2ρ2*tab + (1-2ρ2)*tbw
    tb3 = 2ρ2*(1/me12 + NA) + (1-2ρ2)*tbw
    y, tb1, tb2, tb3, tb
end
plot(Ts[1][1:2], ribbon=(Ts[1][3:4]...,), size=(800,200), color=:gray)
plot!(getindex.(yy, Ref([1,2])), yscale=:log10, lw=2)
plot!(getindex.(yy, Ref([1,3])), yscale=:log10, lw=2)
plot!(getindex.(yy, Ref([1,4])), yscale=:log10, lw=2)
#plot!(getindex.(yy, Ref([1,5])), yscale=:log10, lw=2)
# difference between predictions is minute

# Fst
tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
yy = map(range(0, C, 500)) do x
    y = G*x/C
    me21, me12 = Barriers.me(PP, x)
    _, ta, tb, tab = bidir_ct(NA, NB, me21, me12)
    tbw = tw(NA, NB, me12)
    tb = 2ρ2*tab + (1-2ρ2)*tbw
    y, 1 - tb/tab
end
plot(Ts[4][1:2], ribbon=(Ts[4][3:4]...,), size=(800,200), color=:gray)
plot!(getindex.(yy, Ref([1,2])), lw=2)

# Fst
tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
yy = map(range(0, C, 500)) do x
    y = G*x/C
    me21, me12 = Barriers.me(PP, x)
    tab = 1/me12 + NA
    tbw = tw(NA, NB, me12)
    tb = 2ρ2*tab + (1-2ρ2)*tbw
    y, 1 - tb/tab
end
plot(Ts[4][1:2], ribbon=(Ts[4][3:4]...,), size=(800,200), color=:gray)
plot!(getindex.(yy, Ref([1,2])), lw=2)


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

# Phase-type, no need for it...
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
    Fst, ta, tb, tab
end

