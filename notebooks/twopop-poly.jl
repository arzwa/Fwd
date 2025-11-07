@everywhere using Fwd, Random, StatsBase, ProgressMeter
using Fwd, Random, StatsBase
using Barriers
using Plots, ProgressMeter, WrightDistribution

tab(NA, m) = 1/m + NA
tb(NA, NB, m) = NB*(3 + 2m*NA)/(1 + 2m*NB)
function fst(NA, NB, m) 
    pb = tab(NA, m)
    pw = NA + tb(NA, NB, m)
    (pb - pw)/(pb + pw)
end

function locuseffect(l::Union{HaploidBiLocus,PoissonLocus}, r) 
    l.s == 0. && return 0.0
    l.u/(l.s*(1 + r*(1-l.s)/l.s)^2) 
end
function hkbgs(A::Architecture, R, x)
    logB = 0.0
    for i=1:length(A)
        r = Fwd.recrate(R, x, A.xs[i])
        logB -= locuseffect(A[i], r)
    end
    exp(logB)
end

NA = 500
NB = 500
Ls = 0.25
L1 = 150
s  = Ls/L1
m  = 0.1*s
u  = s/200
L2 = 150
Us = 0.1
sd = 0.01
ud = Us*sd
C  = 1.0   # 1M
G  = 10^8  # 100Mb  -> 1cM/Mb

rng = Random.seed!(28)
dfe = Exponential(s)
ss = rand(rng, dfe, L1)
L = L1 + L2
ys = cumsum(rand(rng, Dirichlet(L+1,10.0)))[1:end-1] .* C
xs = ceil.(Int64, ys .* (G/C))

o = randperm(rng, L)
bgsloci = [PoissonLocus(sd,ud) for _=1:L2]
Aloci = [Fwd.HaploidBiLocus(0.0, 0.0) for _=1:L1]
Bloci = [Fwd.HaploidBiLocus(-ss[i], u) for i=1:L1]
AA = Architecture([Aloci ; bgsloci][o], xs)
AB = Architecture([Bloci ; bgsloci][o], xs)
R  = Fwd.LinearPhysMap(maplength=C, physlength=G)
idx = filter(i-> typeof(AB[i]) <: HaploidBiLocus, 1:length(AB))

xA = [ [ ones(Int, L1) ; rand(rng, Poisson(U/sd), L2)][o] for _=1:NA]
xB = [ [zeros(Int, L1) ; rand(rng, Poisson(U/sd), L2)][o] for _=1:NB]

xx = 1:10000:G
B = map(x->hkbgs(AA, R, x), xx)
plot(xx, B)

nA = collect(1:NA)
nB = collect(1:NB) .+ NA
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 10^5

model, AM = let
    loci = [Barriers.DiploidLocus(-2AB[i].s, 0.5, AB[i].u) for i in idx]
    _ys = ys[idx]
    R = [Fwd.recrate(abs(_ys[i] - _ys[j])) for i=1:L1, j=1:L1]
    A = Barriers.Architecture(loci, _ys, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M)
    AM = AeschbacherModel(m, [AB.loci[i].s for i in idx], _ys)
    EM, AM
end

mx = 0:0.002:C
plot(mx, x->tab(NA, Barriers.me(model,x)), 
    color=:black, ylim=(1000,10ngen), yscale=:log10)
hline!([ngen], ls=:dash)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:salmon)
plot!(twinx(), mx, map(x->fst(NA, NB, Barriers.me(model,x)), mx), 
    color=:gray, ylim=(0,1), framestyle=:default)

# Simulation
pop, ts, qs = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop, G)
    qs=Matrix{Float64}(undef, ngen, L)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        qs[i,:] .= mean(pop.popB.x)
        if i % 50 == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    pop, ts, qs
end

pts = to_tskit(ts)
pth = mkpath("data/2025-09-11/")
pts.dump(joinpath(pth, "ts.ts"))
serialize(joinpath(pth, "qs.jls"), (xs, ss, qs))

# Visualization
p̄ = 1 .- vec(mean(qs, dims=1))[idx]
scatter(model.Ep, p̄, color=:black, ms=2)
plot!(x->x, xlim=(0,1), ylim=(0,1), size=(220,210), 
    xlabel="\$\\mathbb{E}[p]\$", ylabel="\$\\bar{p}\$")

xx, pa, pb, dab = Fwd.diffdiv(ts)
mes = map(1:length(xx)) do i
    a = i == 1 ? 0.0 : xx[i-1]
    b = xx[i]
    a *= C/G
    b *= C/G
    me1 = quadgk(x->Barriers.me(model, x), a, b)[1] ./ (b-a) 
    me2 = quadgk(x->m*Barriers.gff(AM, x), a, b)[1] ./ (b-a) 
    me1, me2
end
me1 = first.(mes)
me2 = last.(mes)
δ = [xx[1]; [xx[i]-xx[i-1] for i=2:length(xx)]]

Bs = map(x->hkbgs(AA, R, x), xx)
NAs = NA .* Bs
NBs = NB .* Bs
plot(xx, pa, line=:steppre, color=:gray, alpha=0.4, yscale=:log10)
plot!(xx, NAs, color=2)
hline!([NA])

p1 = plot(xx, dab, line=:steppre, color=:gray, 
    alpha=0.4, yscale=:log10, ylabel="\$T_{AB}\$")
tab1 = map(i->tab(NAs[i], me1[i]), 1:length(me1))
tab2 = map(i->tab(NAs[i], me2[i]), 1:length(me2))
plot!(xx, tab1)
plot!(xx, tab2)
hline!([ngen], ls=:dash, color=:gray, alpha=0.5)
hline!([sum(dab .* δ) / G])
ix = findall(x->x>=ngen, dab)
scatter!(xx[ix], dab[ix], color=:red, ms=2)

p2 = plot(xx, pb , line=:steppre, color=:gray, alpha=0.4, ylabel="\$T_{B}\$")
hline!([sum(pb .* δ) / G], color=:black, ls=:dash)
tb1 = map(i->tb(NAs[i], NBs[i], me1[i]), 1:length(me1))
tb2 = map(i->tb(NAs[i], NBs[i], me2[i]), 1:length(me2))
plot!(xx, tb1)
plot!(xx, tb2)
plot!(size=(900,200), yscale=:log10, xlabel="map position")
plot(p1, p2, layout=(2,1), size=(900,400), margin=3Plots.mm)

plot(xx, tab1 .- dab)
plot!(xx, tab2 .- dab, ylim=(-1e5,1e5), size=(900,200))

xx, pa, pb, dab = Fwd.diffdiv(ts)
p1 = plot(xx, dab, line=:steppre, color=:gray)
hline!([ngen])
plot!(xx, tab1) 
plot!(xx, tab2)
plot!(size=(900,200), ylim=(0,ngen*1.2))

# coarse model?
Δ  = 0.02
Δs = fill(Δ, ceil(Int64, C/Δ))
zs = [0 ; cumsum(Δs)]
sm = map(1:length(zs)-1) do k
    jdx = filter(i->zs[k] < ys[i] <= zs[k+1], idx)
    length(jdx) == 0 ? 0.0 : mean([-AB[j].s for j in jdx])
end
Xc = fit(Histogram, ys, zs).weights
CM = Barriers.CoarseModel(X=Xc, Δ=Δs, s=sm, m=m, u=u, λ=1/NB)
mec = Barriers.gff(CM) .* m
p1 = plot(0:Δ:C, [mec[1] ; mec], line=:steppre, 
    title="\$m_e\$", color=:teal, lw=1.5, xlabel="map position")
plot!(mx, x->Barriers.me(model,x), color=:black, alpha=0.5)
p2 = plot(0:Δ:C, [Xc[1]; Xc], line=:steppre, color=:firebrick, fill=true,
    fillalpha=0.3, title="number of selected sites \$X_i\$")
plot(p2, p1, layout=(2,1))

xx, pa, pb, dab = Fwd.diffdiv(ts)
plot(xx ./ G, dab, line=:steppre, color=:lightgray, label="simulation")
plot!(mx, x->tab(NA, Barriers.me(model, x)), 
    color=:black, yscale=:log10, label="approx. 1 (diffusion)")
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, 
    color=:teal, label="approx. 2 (Aeschbacher et al.)")
plot!(0:Δ:C, tab.(NA, [mec[1];mec]), 
    color=:teal, line=:steppre, label="coarse model")
plot!(ylabel="\$T_{AB}\$", xlabel="\$x\$", legend=:topleft,
    size=(700,220), margin=4Plots.mm)
sticks!(twinx(), ys[idx], ss, ylabel="\$s\$", framestyle=:default, 
    color=:firebrick)

plot(x->m*Barriers.gff(AM, x), mx, 
    color=:teal, label="approx. 2 (Aeschbacher et al.)")
plot!(mx, x->Barriers.me(model,x), color=:black, )


# simulate data
pts = Fwd.tskit.load("data/2025-09-09.ts")
using PyCall
msprime = pyimport("msprime")
ts = msprime.sim_mutations(pts, rate=0.1, 
    random_seed=22, discrete_genome=false)

G = ts.genotype_matrix()


# --- reps ---------------------------------------------------------------
let mpop=mpop
    res = pmap(1:10) do rep
        seed = rand(1:2^32)
        rng = Random.seed!(seed)
        pop = deepcopy(mpop)
        ts = Fwd.init_ts(pop, G) 
        qs = Matrix{Float64}(undef, ngen, L)
        for i=1:ngen
            pop = Fwd.generation!(rng, pop, ts);
            qs[i,:] .= mean(pop.popB.x)
            if i % 20 == 0 
                pop, ts = Fwd.simplify!(pop, ts)
            end
        end
        seed, pop, ts, qs
    end
end

dv = map(res) do (seed, pop, ts)
    xx, pa, pb, dab = Fwd.diffdiv(ts)
    xx[2:end], pb, dab 
end
q̄ = vec(mean(vcat(map(Q->Q[10NB:end,:], last.(res))...), dims=1))

# simulation
x, y = Fwd.summarize_wins(first.(dv), getindex.(dv, 3))

p1 = plot(x, vec(mean(y, dims=1)), label="simulations (\$n=10\$)",
    line=:steppre, yscale=:log10, color=:lightgray) 
# diffusion model
mx = 0:0.002:C
plot!(x->tab(NA, Barriers.me(model, x)), mx, color=:black, lw=1, label="diffusion")
# Petry (closest barrier)
_gff(r, s) = r/(r+s)
tt = map(x->tab(NA, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
#plot!(mx, tt, color=:orange, lw=1, label="Petry (closest barrier)")
# Aeschbacher
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:firebrick, label="Aeschbacher")
plot!(0:Δ:C, tab.(NA, [mec[1];mec]), color=:teal, line=:steppre, label="coarse model")
plot!(title="\$T_{AB}, Ls=$Ls, L=$L, N_e s=$(NB*s), m/s=$(m/s)\$", 
    legend=:outertopright, xlabel="map position (\$M\$)", 
    size=(800,250))

p2 = scatter(1 .- q̄, model.Ep, color=:black, ms=2, 
    xlabel="simulation", ylabel="prediction", title="freq. locally beneficial allele")
plot!(x->x, color=:black, size=(300,300), xlim=(0,1), ylim=(0,1))

plot(p1, p2, size=(1100,280), layout=grid(1,2,widths=[0.78,0.22]), margin=5Plots.mm)



x, y = Fwd.summarize_wins(first.(dv), getindex.(dv, 2))
p2 = plot(x, vec(mean(y, dims=1)), line=:steppre, yscale=:log10, color=:lightgray) 
mx = 0:0.002:C
tt = map(x->tb(NA, NB, Barriers.me(model, x)), mx)
plot!(mx, tt, color=:black, lw=2)
tt = map(x->tb(NA, NB, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
plot!(mx, tt, color=:orange, lw=2)
plot!(x->tb(NA, NB, m*Barriers.gff(AM, x)), mx, color=:cyan)
title!("\$T_{B}\$")
plot(p1, p2, size=(700,250), xlabel="map position (M)", margin=3Plots.mm)


