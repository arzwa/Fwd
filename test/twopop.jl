@everywhere using Fwd, Random, StatsBase
using Plots, ProgressMeter, WrightDistribution

NA = 500
NB = 500
L = 2
C = 1.0
s = 0.05
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
ngen = 100*(NB+NA)
res = pmap(1:10) do rep
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
plot!(mx, tt, color=:black, lw=2)
_gff(r, s) = r/(r+s)
tt = map(x->tab(NA, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
plot!(mx, tt, color=:orange, lw=2)
AM = AeschbacherModel(m, fill(-s, L), xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:cyan)
plot!(title="\$T_{AB}\$")
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

