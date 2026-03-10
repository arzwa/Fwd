using Fwd, Distributions, Random, Plots, Parameters
using Barriers
import TreeSequences as TS
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
NA   = 1
Ns   = 5. 
NB   = ceil(Int64, Ns/s̄)
loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
mBC = 0.
mCB = 1.5s̄
NB  = ceil(Int64, Ns/s̄)
NA  = NB
NC  = 2NB
C   = 0.5
G   = C*100*10^6   # (1cM/Mb)
xs  = ceil.(Int64, ys .* G)
R   = LinearPhysMap(C=C, G=G)

function secondary_contact_model(
        xs, ss, R, u, mBC, mCB, NA, NB, NC, ngen; kwargs...)
    # Ancestral population: this one starts from `00...0`.
    # the ancestral population evolves neutrally at all L loci
    L = length(xs)
    arch = Architecture([BiAllelic(u) for _=1:L], xs, R)
    gpmA = GPMap([HaploidLocus(0.0, i) for i=1:L])
    nA = collect(1:NA)
    xA = [rand(Bool, L) for _=1:NA]
    popA = WFPopulation(ploidy=Haploid(), 
        N=NA, arch=arch, gpm=gpmA, x=xA, nodes=nA)
    popA, ts, qs1 = Fwd.simulate!(popA, init_ts(popA), ngen[1], 
        pop->mean(pop.x); kwargs...)
    pops = Fwd.splitpop(popA, ts, (NB, NC))
    gpmB = GPMap([HaploidLocus( ss[i], i) for i=1:L])
    gpmC = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
    pops[1] = reconstruct(pops[1], gpm=gpmB)
    pops[2] = reconstruct(pops[2], gpm=gpmC)
    # Allopatric phase
    mpop = Fwd.MetaPop(pops, [0. 0.; 0. 0.])
    mpop, ts, qs2 = Fwd.simulate!(mpop, ts, ngen[2], 
        pop->(mean(pop.P[1].x), mean(pop.P[2].x)); kwargs...)
    # Secondary contact phase
    mpop.M[1,2] = mBC
    mpop.M[2,1] = mCB
    mpop, ts, qs3 = Fwd.simulate!(mpop, ts, ngen[3], 
        pop->(mean(pop.P[1].x), mean(pop.P[2].x)); kwargs...)
    PA = pcat(qs1...) 
    PB2 = pcat(first.(qs2)...)
    PB3 = pcat(first.(qs3)...)
    PC2 = pcat(last.( qs2)...)
    PC3 = pcat(last.( qs3)...)
    Ps =  (PA, PB2, PB3, PC2, PC3)
    return mpop, ts, Ps
end

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

ngen = [10*max(NB,NC) for _=1:3]
nrep = 10
res = map(1:nrep) do k
    @info k
    secondary_contact_model(xs, ss, R, u, mBC, mCB, NA, NB, NC, ngen; 
        simplify=100, every=ngen[1]÷200)
end

# Allele frequency divergence
div = map(last.(res)) do (PA, PB1, PB2, PC1, PC2)
    pb = vec(mean(PB2, dims=1))
    pc = vec(mean(PC2, dims=1))
    pb .- pc
end 

tss = getindex.(res, 2)
xx, tb = estimate_coaltimes(tss, idx=2)
xx, tc = estimate_coaltimes(tss, idx=3)
xx, tbc = estimate_coaltimes(tss, idx=4)

fst = 1 .- ((tb .+ tc) ./ 2) ./ tbc
plot(xx, fst, size=(800,200), color=:gray, alpha=0.5, lw=0.5)
sticks!(xs, mean(div), )

function pgfst(NA, NB, NC, t1, t2, mCB, mBC)
    # t1 is time of SC (back in time)
    # t2 is time of split (back in time)
    N = NB
    demography = pg.Demography(
        pop_sizes=Dict("B"=>NB/N, "C"=>NC/N),
        migration_rates=Dict(
            ("B","C")=>mCB*N, 
            ("C","B")=>mBC*N))
    demography.add_event(pg.PopulationSplit(
        derived="C", ancestral="B", time=t2/N))
    demography.add_event(pg.PopSizeChange(
        pop="B", time=t2/N, size=NA/N))
    demography.add_event(pg.MigrationRateChange(
        source="B", dest="C", time=t1/N, rate=0.))
    demography.add_event(pg.MigrationRateChange(
        source="C", dest="B", time=t1/N, rate=0.))
    coal = pg.Coalescent(
        n = Dict("B"=>1, "C"=>1),
        demography = demography)
    tbc = coal.tree_height.mean * N
    coal = pg.Coalescent(
        n = Dict("B"=>2, "C"=>0),
        demography = demography)
    tb = coal.tree_height.mean * N
    coal = pg.Coalescent(
        n = Dict("C"=>2, "B"=>0),
        demography = demography)
    tc = coal.tree_height.mean * N
    Fst = 1 - ((tb + tc) / 2) / tbc
    Fst, tbc, tb, tc
end
    
BP = Equilibrium(BPModel(m=mCB, xs=ys*C, Ne=NB, u=u, s=ss))
mes = map(x->(G*x/C, Barriers.me(BP,x)), range(0,C,500))
yy = map(mes) do (x, me)
    x, pgfst(NA, NB, NC, ngen[3], ngen[3]+ngen[2], me, mBC)
end

# Fst scan plot
plot(xx, fst, ylim=(0,1), color=:gray, alpha=0.8, lw=0.5, size=(800,200))
sticks!(xs, mean(div), lw=3, color=:firebrick, alpha=0.3)
plot!(first.(yy), first.(last.(yy)), lw=1, color=:black)
plot!(xlabel="map position", ylabel="\$F_\\mathrm{ST}\$", margin=5Plots.mm)

PP = map(res) do (_,_,Ps)
    PA, PB1, PB2, PC1, PC2 = Ps
    vcat(PA, PB1, PB2), vcat(PA, PC1, PC2)
end

# Example allele frequency trajectories
map(enumerate(5:5:L)) do (j,k)
    pl = plot(title=@sprintf("locus %d (\$s = %.1e\$)", k, ss[k]))
    plot!(ylabel = j % 2 == 1 ? "\$p\$" : "")
    plot!(xlabel = j >= 9  ? "\$t\$" : "")
    tt = 1:100:100*size(PP[1][1], 1)
    map(PP) do P
        plot!(tt, P[1][:,k], color=1, alpha=0.5)
        plot!(tt, P[2][:,k], color=2, alpha=0.5)
    end
    hline!([mean(mean(first.(PP))[end-150:end,k])], color=:black, ls=:dot)
    hline!([BP.Ep[k]], color=:black, ls=:dash)
    xts = cumsum([0; ngen])
    plot!(rectangle(xts[2], xts[3], 0, 1), lw=0, color=:gray, alpha=0.2, fill=true)
    plot!(xlim=(0,Inf), xticks=xts, ylim=(0,1))
end |> x->plot(x..., layout=(5,2), size=(500,850), 
    right_margin=5Plots.mm, left_margin=3Plots.mm)

divest = map(PP) do (PB,PC)
    vec(mean(PB[end-100:end,:], dims=1)) .- vec(mean(PC[end-100:end,:], dims=1))
end |> mean

scatter(divest, BP.Ep)


# -------------
map(1:20) do k
    plot(PB[:,k])
    plot!(PC[:,k])
end |> x->plot(x..., layout=(5,4), size=(700,600))

pb = mean(PB, dims=1) |> vec
pc = mean(PC, dims=1) |> vec
scatter(pb, pc, ylim=(0,1), xlim=(0,1), 
    size=(300,280), ms=2, color=:black,
    title="\$N_Bm_{CB} = $(NB*mCB), N_Cm_{BC} = $(NC*mBC)\$",
    xlabel="\$p_B\$", ylabel="\$p_C\$")

k = 8
stephist( PB[:,k], bins=0:0.05:1.05)
stephist!(PC[:,k], bins=0:0.05:1.05)


plot(xx, tbc, size=(1000,200), )
plot!(xx, tb)
plot!(xx, tc)
vline!(xs)

# Fst
fst = 1 .- ((tb .+ tc) ./ 2) ./ tbc
P1 = plot(xx, fst, size=(1000,200), color=:black, framestyle=:default, ylabel="\$F_{ST}\$", xlabel="map position", margin=5Plots.mm, alpha=0.2)
sticks!(twinx(), xs, pb .- pc, lw=5, color=:firebrick, alpha=0.3, 
    ylim=(0,1), framestyle=:default)


P3 = scatter(pb, BP.Ep, ms=2, color=:black, xlabel="simulation", ylabel="prediction")
plot!(x->x, xlim=(0,1), ylim=(0,1))

# should do proper SC prediction but with mₑ 
function TB(N1, N2, m12, ng2, ng3)
    M = N2*m12
    x = N1/N2
    τ = ng3/2N2
    τa = (ng3 + ng2)/2N2
    1 + 2/M + exp(-τ*M/2)*(τa - 2/M) + (1-x)*(
        (M/(M-2)) * exp(-τ*M/2 - τa) - 
        (M/(M-2)) * exp(-τ - τa) - exp(-τ*M/2))
end

TB(NC, NB, mCB, ngen2, ngen3)

mes = map(x->(G*x/C, 1/(TB(NC, NB, Barriers.me(BP,x), ngen2, ngen3)*2NB)), range(0,C,500))
P2 = plot(xx, 1 ./ tbc, size=(1000,200), color=:gray, alpha=0.6)
plot!(mes, yscale=:log10, lw=2, color=:black, xlabel="map position",ylabel="\$1/T_{b}\$", margin=5Plots.mm)
sticks!(twinx(), xs, pb .- pc, lw=5, color=:firebrick, alpha=0.3, 
    ylim=(0,1), framestyle=:default)

P4 = plot(P1, P2, layout=(2,1), size=(1100,300))
plot(P4,P3,layout=grid(1,2,widths=[0.75,0.25]))


# Phase gen?

fst = 1 .- ((tb .+ tc) ./ 2) ./ tbc
P1 = plot(xx, fst, size=(1000,200), color=:black, framestyle=:default, ylabel="\$F_{ST}\$", xlabel="map position", margin=5Plots.mm, alpha=0.2)
sticks!(twinx(), xs, pb .- pc, lw=5, color=:firebrick, alpha=0.3, 
    ylim=(0,1), framestyle=:default)
plot!(first.(yy), first.(last.(yy)), lw=2, color=:firebrick)


# Bidirectional mₑ prediction does not seem to work...
# B is the population in which the alleles are favored
# C is the population in which the alleles are deleterious (p ~ 0)
BP2 = Barriers.BPTwoPop(
    m12 = mCB,
    m21 = mBC,
    s1  = ss,
    s2  = ss,
    xs  = ys * C,
    Ne1 = 1.0NC,
    Ne2 = 1.0NB,
    u   = u)

PP = Equilibrium(BP2, α=0.05)

scatter(BP.Ep, PP[end,L+1:end])
plot!(x->x)

map(1:L) do k
    plot(PP[:,k])
    plot!(PP[:,L+k])
end |> x->plot(x..., size=(900,700))

P1 = scatter(pb, PP[end,L+1:2L])
plot!(x->x)
P2 = scatter(pc, PP[end,1:L])
plot!(x->x)
plot(P1, P2)
