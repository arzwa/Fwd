@everywhere using Fwd, Random, StatsBase
using Fwd, Random, StatsBase
using Barriers
using Plots, ProgressMeter, WrightDistribution
tab(NA, m) = 1/m + NA
tb(NA, NB, m) = NB*(3 + 2m*NA)/(1 + 2m*NB)

NA = 500
NB = 500
Ls = 0.4
L  = 35
s  = Ls/L
C  = 1.0
m  = 0.5*s
@info NB*s
u  = s/200
rng = Random.seed!(282)
dfe = Exponential(s)
xs = cumsum(rand(rng, Dirichlet(L+1,10.0)))[1:end-1] .* C
AA = Architecture([Fwd.HaploidBiLocus(0.0, 0.0) for _=1:L], xs)
AB = Architecture([Fwd.HaploidBiLocus(-rand(rng, dfe), u) for _=1:L], xs)
ss = [-AB.loci[i].s for i=1:L]
R  = LinearMap(C)
xA = [ ones(Bool, L) for _=1:NA]
xB = [zeros(Bool, L) for _=1:NB]
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 5*10^5
pop = deepcopy(mpop)
ts = Fwd.init_ts(pop, C) 
model = let
    recmap = Barriers.linearmap(100, C)
    loci = [Barriers.DiploidLocus(2ss[i], 0.5, u) for i=1:L]
    R = [Fwd.recrate(abs(xs[i] - xs[j])) for i=1:L, j=1:L]
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M)
end
mx = 0:0.002:C
plot(mx, x->tab(NA, Barriers.me(model,x)), 
    color=:black, ylim=(1000,10ngen), yscale=:log10)
hline!([ngen], ls=:dash)
AM = AeschbacherModel(m, [AB.loci[i].s for i=1:L], xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:salmon)
plot!(twinx(), mx, x->Barriers.me(model,x), 
    color=:cyan, framestyle=:default)

@showprogress for i=1:ngen
    pop = Fwd.generation!(rng, pop, ts);
    if i % 10 == 0 
        pop, ts = Fwd.simplify!(pop, ts)
    end
end

xx, pa, pb, dab = Fwd.diffdiv(ts)
plot(xx, dab, line=:steppre, color=:lightgray)
plot!(mx, x->tab(NA, Barriers.me(model,x)), color=:black, )
AM = AeschbacherModel(m, [AB.loci[i].s for i=1:L], xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:salmon)
plot!(size=(900,200), yscale=:log10)

scatter(mean(pop.popB.x), 1 .- model.Ep, color=:black, ms=2)
plot!(x->x, color=:gray, size=(210,200))

# coarse model?
Δ  = 0.005
Δs = fill(Δ, ceil(Int64, C/Δ))
zs = [0 ; cumsum(Δs)]
sm = map(1:length(zs)-1) do k
    idx = filter(i->zs[k]< xs[i] <= zs[k+1], 1:length(xs))
    length(idx) == 0 ? 0.0 : mean(ss[idx])
end
Xc = fit(Histogram, xs, zs).weights
CM = Barriers.CoarseModel(X=Xc, Δ=Δs, s=sm, m=m, u=u, λ=0.0)
mec = Barriers.gff(CM) .* m
plot(0:Δ:C, [mec[1] ; mec], line=:steppre)
plot!(mx, x->Barriers.me(model,x), color=:black, )

xx, pa, pb, dab = Fwd.diffdiv(ts)
plot(xx, dab, line=:steppre, color=:lightgray, label="simulation")
plot!(mx, x->tab(NA, Barriers.me(model,x)), 
    color=:black, yscale=:log10, label="approx. 1 (diffusion)")
AM = AeschbacherModel(m, [AB.loci[i].s for i=1:L], xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, 
    color=:teal, label="approx. 2 (Aeschbacher et al.)")
#plot!(0:Δ:C, tab.(NA, [mec[1];mec]), 
#    color=:teal, line=:steppre, label="coarse model")
plot!(ylabel="\$T_{AB}\$", xlabel="\$x\$", legend=:topleft,
    size=(700,220), margin=4Plots.mm)
sticks!(twinx(), xs, ss, ylabel="\$s\$", framestyle=:default, 
    color=:firebrick)

plot(x->m*Barriers.gff(AM, x), mx, 
    color=:teal, label="approx. 2 (Aeschbacher et al.)")
plot!(mx, x->Barriers.me(model,x), color=:black, )

pts = to_tskit(ts)
pts.dump("data/2025-09-09.ts")

# simulate data
pts = Fwd.tskit.load("data/2025-09-09.ts")
using PyCall
msprime = pyimport("msprime")
ts = msprime.sim_mutations(pts, rate=0.1, 
    random_seed=22, discrete_genome=false)

G = ts.genotype_matrix()


# --- reps ---------------------------------------------------------------
res = @showprogress pmap(1:10) do rep
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    pop = deepcopy(mpop)
    ts = Fwd.init_ts(pop, C) 
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


# ---------------------------------------------------------

NA = 100
NB = 100
L = 25
C = 1.0
s = 0.02
m = 0.001
u = s/200
xs = collect(C/2L:C/L:C)
AA = Architecture([Fwd.HaploidBiLocus(0.0, 0.0) for _=1:L], xs)
AB = Architecture([Fwd.HaploidBiLocus( -s, u  ) for _=1:L], xs)
R  = LinearMap(C)
xA = [ ones(Bool, L) for _=1:NA]
xB = [zeros(Bool, L) for _=1:NB]
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
rng = Random.seed!(15)
ngen = 1000*(NB+NA)
res = pmap(1:24) do rep
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    pop = deepcopy(mpop)
    ts = Fwd.init_ts(pop, C) 
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        if i % 100 == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    seed, pop, ts
end

dv = map(res) do (seed, pop, ts)
    xx, pa, pb, dab = Fwd.diffdiv(ts)
    xx[2:end], pb, dab 
end

using Barriers
model = let
    recmap = Barriers.linearmap(100, C)
    loci = fill(Barriers.DiploidLocus(2s, 0.5, s/1000), L)
    R = [Fwd.recrate(abs(xs[i] - xs[j])) for i=1:L, j=1:L]
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M)
end

tab(NA, m) = 1/m + NA
tb(NA, NB, m) = (3NB − 4NB*m + 2NA*NB*m + m^2*NB − m^2*NA*NB)/(1 − 2m + 2NB*m + m^2 − m^2*NB)

x, y = Fwd.summarize_wins(first.(dv), getindex.(dv, 3))
p1 = plot(x, vec(mean(y, dims=1)), line=:steppre, yscale=:log10, color=:lightgray) 
mx = 0:0.002:C
tt = map(x->tab(NA, Barriers.me(model, x)), mx)
plot!(mx, tt, color=:black)
_gff(r, s) = r/(r+s)
tt = map(x->tab(NA, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
plot!(mx, tt, color=:orange)
AM = AeschbacherModel(m, fill(-s, L), xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:cyan)
plot!(title="\$T_{AB}\$")
x, y = Fwd.summarize_wins(first.(dv), getindex.(dv, 2))
p2 = plot(x, vec(mean(y, dims=1)), line=:steppre, yscale=:log10, color=:lightgray) 
mx = 0:0.002:C
tt = map(x->tb(NA, NB, Barriers.me(model, x)), mx)
plot!(mx, tt, color=:black)
tt = map(x->tb(NA, NB, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
plot!(mx, tt, color=:orange)
plot!(x->tb(NA, NB, m*Barriers.gff(AM, x)), mx, color=:cyan)
title!("\$T_{B}\$")
plot(p1, p2, size=(700,250), xlabel="map position (M)", margin=3Plots.mm)

