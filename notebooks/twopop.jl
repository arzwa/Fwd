@everywhere using Fwd, Random, StatsBase
using Fwd, Random, StatsBase, Parameters
using ProgressMeter, WrightDistribution
using Plots

NA = 1
NB = 500
L = 1
C = 1.0
s = 0.05
m = 0.001
u = s/500
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
ngen = 50NB

pop = deepcopy(mpop)
ts = Fwd.init_ts(pop, C) 
@showprogress for i=1:ngen
    pop = Fwd.generation!(rng, pop, ts);
    if i % 100 == 0 
        pop, ts = Fwd.simplify!(pop, ts)
    end
end

# ---------------------------------------------------------------
# try to get an idea of block lengths with migrant ancestry
# generated unsimplified ts for additional `n` gens
function trace_descendants(ts, node, C, tmax)
    segments = Tuple{Int,Float64,Float64}[]
    t = time(ts[node])
    function _recurse(n, x0, x1)
        time(ts[n]) - t >= tmax && return
        for ei in ts.children[n]
            edge = ts.edges[ei]
            # only consider sink population descendants
            Fwd.population(ts[edge.child]) != 2 && continue
            y0, y1 = edge.left, edge.rght
            if y0 <= x1 && x0 <= y1  # overlapping with parent segment
                z0 = max(x0,y0)
                z1 = min(x1,y1)
                push!(segments, (time(ts[edge.child]) - t, z0, z1))
                _recurse(edge.child, z0, z1)
            end
        end
    end
    _recurse(node, 0.0, C)
    segments
end

function segdensity(segments, bins)
    cm = countmap(segments)
    xs = map(2:length(bins)) do i
        x0, x1 = bins[i-1], bins[i]
        n = 0
        for (k,v) in cm
            _, y0, y1 = k
            if y0 <= x1 && x0 <= y1
                n += v
            end
        end
        n
    end
end

rng = Random.seed!(17)
nrep = 10_000
tps = [1,2,5,10,20,50,100,200]
res = @showprogress map(1:nrep) do j
    dpop = reconstruct(deepcopy(pop), m=0.0)
    dts = deepcopy(ts)
    k = rand(rng, 1:dpop.popB.N)
    nk = Fwd.migrate!(dpop.popA, dpop.popB, 1, k)
    W = Fwd.eval_fitness(dpop.popB)
    for i=1:tps[end]
        dpop = Fwd.generation!(rng, dpop, dts);
    end
    map(tps) do gen
        segments = trace_descendants(dts, nk, C, gen) |> sort
        segs = filter(x->x[1] == gen, segments)
    end
end

xx = 0:0.005:C
sds = map(res) do rep
    map(X->segdensity(X, xx), rep)
end

P1=plot()
ps = map(zip(tps, sum(sds))) do (tp, X)
    y = X ./ nrep
    y_ = y ./ step(xx)
    plot!(P1, xx, [y_[1];y_], line=:steppre, label="\$t=$tp\$", ylabel="# segments")
    plot(xx, [y[1];y], line=:steppre, title="\$t=$tp\$", ylabel="# segments")
end
plot(ps..., margin=4Plots.mm, layout=(4,2), size=(700,900))

plot(P1, ylabel="#segments/M", yscale=:log10, legend=:bottomright)

idx =findall(x->!isempty(x[1]), res)

j = idx[4291]
map(enumerate(res[j])) do (i,y)
    p = plot(xlim=(0,C), framestyle=:box, yticks=false, 
        ylim=(0,maximum(map(length, y), init=0)), title="\$t=$(tps[i])\$")
    a = 0.05
    for (i,x) in enumerate(y)
        plot!([x[2], x[3]], a*[i,i], color=:black, lw=1.5)
    end
    p
end |> x->plot(x..., layout=(4,2), size=(400,600))


# ---------------------------------------------------------------

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

#  TAB
#  simulation
x, y = Fwd.summarize_wins(first.(dv), getindex.(dv, 3))
p1 = plot(x, vec(mean(y, dims=1)), line=:steppre, yscale=:log10, color=:lightgray) 
#  polygenic me prediction
mx = 0:0.002:C
tt = map(x->tab(NA, Barriers.me(model, x)), mx)
plot!(mx, tt, color=:black, lw=2)
#  Petry me
_gff(r, s) = r/(r+s)
tt = map(x->tab(NA, m*_gff(Fwd.recrate(minimum(abs.(x .- xs))), s)), mx)
plot!(mx, tt, color=:orange, lw=2)
#  Aeschbacher
AM = AeschbacherModel(m, fill(-s, L), xs)
plot!(x->tab(NA, m*Barriers.gff(AM, x)), mx, color=:cyan)
plot!(title="\$T_{AB}\$")
#  TB
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

