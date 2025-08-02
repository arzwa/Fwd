using Fwd, Random, StatsBase, Plots, ProgressMeter, WrightDistribution

NA = 500
NB = 800
L = 10
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
ts = Fwd.init_ts(mpop, C) 
rng = Random.seed!(15)
ngen = 100*(NB+NA)
@showprogress for i=1:ngen
    mpop = Fwd.generation!(rng, mpop, ts);
    if i % 100 == 0 
        mpop, ts = Fwd.simplify!(mpop, ts)
    end
end
ts1 = simplify(ts, Fwd.leaves(ts, 1))
ts2 = simplify(ts, Fwd.leaves(ts, 2))

map([(ts1,NA), (ts2,NB)]) do (ts, N)
    xts = reconstruct(ts, nodes=[Fwd.Node(time(n), 0) for n in ts.nodes])
    xs, hs = Fwd.theights(xts)
    plot(xs, hs, line=:steppre)
    hline!([mean(filter(!isnan,hs))], color=:black)
    hline!([2N], ls=:dash, color=:salmon)
end |> x->plot(x..., size=(600,200))

tab(NA, m) = 1/m + NA
tb(NA, NB, m) = (3NB − 4NB*m + 2NA*NB*m + m^2*NB − m^2*NA*NB)/(1 − 2m + 2NB*m + m^2 − m^2*NB)
xx, pa, pb, dab = Fwd.diffdiv(ts)
p1 = plot( xx[2:end], pa)
hline!([mean(pa)], color=:black)
hline!(xx[2:end], [NA], color=:salmon, ls=:dash, lw=2)
p2 = plot(xx[2:end], pb)
hline!([mean(pb)], color=:black)
hline!([tb(NA, NB, m*exp(-2L*s))], ls=:dash, lw=2, color=:salmon)
hline!([tb(NA, NB, m)], ls=:dash, lw=2, color=:cyan)
vline!(xs, color=:black, ls=:dash)
p3 = plot(xx[2:end], dab, yscale=:log10)
hline!([mean(dab)], color=:black)
hline!([tab(NA, m)], color=:salmon, ls=:dash, lw=2)
hline!([ngen])
vline!(xs,color=:black,ls=:dash)
plot(p1, p2, p3, size=(800,200), layout=(1,3))

plot(xx[2:end], dab, size=(900,200), yscale=:log10)
hline!([ngen, tab(NA, m)])

using Barriers
model = let
    recmap = Barriers.linearmap(100, C)
    loci = fill(Barriers.DiploidLocus(2s, 0.5, s/1000), L)
    R = [Fwd.recrate(abs(xs[i] - xs[j])) for i=1:L, j=1:L]
    A = Barriers.Architecture(loci, xs, R)
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=NB, mode=1)
    EM = Barriers.Equilibrium(M)
end

mx = 0:0.01:C
plot(mx, x->Barriers.me(model, x))
tt = map(x->tab(NA, Barriers.me(model, x)), mx)
plot(xx[2:end], dab, size=(900,200))
plot!(mx, tt, color=:black, lw=2)
hline!([ngen, tab(NA, m)])

tt = map(x->tb(NA, NB, Barriers.me(model, x)), mx)
plot(xx[2:end], pb, size=(900,200))
plot!(mx, tt, color=:black, lw=2)


#----
pts = to_tskit(rts)
#pts.simplify(pts.samples(), keep_input_roots=true)

smpl = [pts.samples(population=0)[1:5]; pts.samples(population=1)[1:5]]
pts.simplify(smpl, keep_input_roots=true).draw_text() |> print
draw_text(simplify(rts, [1:5; 2NA+1:2NA+5], keep_roots=true)) |> print


demography = msprime.Demography.from_tree_sequence(pts)
demography.populations[1].initial_size=NA
demography.populations[2].initial_size=NB
demography.set_migration_rate(1, 0, 0.5)
cts = msprime.sim_ancestry(
    demography=demography,
    initial_state=pts,
    discrete_genome=false,
    ploidy=2,
    random_seed=1)

x, dab, da, db = diffdiv(cts)
p1 = plot(x, dab, line=:steppre)
p2 = plot(x, da , line=:steppre)
p3 = plot(x, db , line=:steppre)
x, dab, da, db = diffdiv(pts)
plot!(p1, x, dab, line=:steppre)
plot!(p2, x, da , line=:steppre)
plot!(p3, x, db , line=:steppre)
plot(p1,p2,p3,layout=(1,3), size=(700,200))

theights(pts)
theights(ppts)




# -----------------------------------------------------------------------
    
NA = 100
NB = 500
L = 1
C = 0.5
s = 0.05
m = 0.005
h = 0.5
u = s*h/200
xs = collect(C/2L:C/L:C)
AA = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
AB = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u  ) for _=1:L], xs)
R  = LinearMap(C)
xA = [ones(Bool, L) for _=1:2NA]
xB = [zeros(Bool, L) for _=1:2NB]
nA = collect(1:2NA)
nB = collect(1:2NB) .+ 2NA
popA = Fwd.DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=deepcopy(xA), nodes=nA)
popB = Fwd.DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=deepcopy(xB), nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ts = Fwd.init_ts(mpop, C) 
qs  = Vector{Float64}[]
#mpop = Fwd.generation!(rng, mpop, ts)
#push!(qs, mean(mpop.popB.x))

rng = Random.seed!(1)
ngen = 100NB
for i=1:ngen
    mpop = Fwd.generation!(rng, mpop, ts)
    push!(qs, mean(mpop.popB.x))
    if i % NB == 0
        @info "$i/$ngen" 
        mpop, ts = Fwd.simplify!(mpop, ts)
    end
end
rts = reverse_relabel(ts)
pts = Fwd.to_tskit(rts)

stephist(hcat(qs...)', norm=true)
d = Wright(-2NB*s, 2NB*u, 2NB*(m + u), h)
plot!(0:0.001:1, p->pdf(d,1-p))


hs = theights(pts)
ix = findall(x->length(x) == 1, hs)
xs = collect(pts.breakpoints())[ix .+ 1]
plot(xs, vcat(hs[ix]...), line=:steppre, size=(700,200), yscale=:log10)
vline!(AB.xs)

demography = msprime.Demography.from_tree_sequence(pts)
demography.populations[1].initial_size=NA
demography.populations[2].initial_size=NB
demography.set_migration_rate(1, 0, 0.5)

cts = msprime.sim_ancestry(
    demography=demography,
    initial_state=pts,
    discrete_genome=false,
    ploidy=1,
    random_seed=7)

from_tskit(cts).nodes

x, dab, da, db = diffdiv(cts)
P1 = plot(x, dab, line=:steppre, title="AB", legend=false, yscale=:log10)
hline!([1/m + 2NA], label="")
vline!([AA.xs], color=:black)
P2 = plot(x, da, line=:steppre, title="A")
hline!([2NA])
P3 = plot(x, db, line=:steppre, title="B")
plot(P1, P2, P3, layout=(1,3), legend=false, size=(900,200))


# -------------------------------
nrep = 20
res = map(1:nrep) do k
    popA = Fwd.DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=deepcopy(xA), nodes=nA)
    popB = Fwd.DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=deepcopy(xB), nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    ts = Fwd.init_ts(mpop, C) 
    @info k
    qs  = Vector{Float64}[]
    for i=1:10NB
        mpop = Fwd.generation!(rng, mpop, ts)
        push!(qs, mean(mpop.popB.x))
        if i % 100 == 0
            #@info i
            mpop, ts = Fwd.simplify!(mpop, ts)
        end
    end
    rts = reverse_relabel(ts)
    pts = Fwd.to_tskit(rts)
    s0 = pts.samples(population=0)
    s1 = pts.samples(population=1)
    xs = collect(pts.breakpoints())
    div = pts.divergence([s0, s1], mode="branch", windows=xs) ./ 2
    piA = pts.diversity(s0, mode="branch", windows=xs) ./ 2
    piB = pts.diversity(s1, mode="branch", windows=xs) ./ 2
    (xs, div, piA, piB, qs, rts) 
end

X, Y = Fwd.summarize_wins(map(y->y[2:end], getindex.(res,1)), getindex.(res, 2))
plot(X, vec(mean(Y,dims=1)))


pts.simplify(90:100).at_index(20).draw_text() |> print

pts.divergence([collect(0:2NA-1), collect(2NA:(2NA + 2NB -1))], 
    windows=collect(pts.breakpoints()))

pts.diversity(collect(0:(2NA + 2NB-1)), 
    windows=collect(pts.breakpoints()))
 
# ----------------
rng = Random.seed!(32)
res = map(1:100) do l
    @info l
    ts = Fwd.init_ts(2NA + 2NB, C) 
    xA = [ones(Bool, L) for _=1:2NA]
    xB = [zeros(Bool, L) for _=1:2NB]
    popA = Fwd.DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=xA, nodes=collect(1:2NA))
    popB = Fwd.DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=xB, nodes=collect(2NA+1:2NA+2NB))
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    q = Vector{Float64}[]
    for i=1:50NB
        mpop = Fwd.generation!(rng, mpop, ts)
        push!(q, mean(mpop.popB.x))
        if i % 100 == 0
            mpop, ts = Fwd.simplify!(mpop, ts)
        end
    end
    ts, q
end

ts = res[1][1]
es = [filter(x->x.child == n, ts.edges) for n in 1000:1100]
ls = map(ex->sum([x.rght-x.left for x in ex]), es)


pyts = Fwd.to_tskit(rts)
xs = collect(pyts.breakpoints())
ss = [collect(1:2NA), collect(1:2NB)]
ps = pyts.divergence(sample_sets=ss, mode="branch", windows=xs)
xs[2:end], ps, q

Xs, Ys = Fwd.summarize_wins(first.(res), getindex.(res, 2))
plot(Xs, vec(mean(Ys, dims=1)))


plot(q); hline!([m/(s*h) mean(q)])
    
tt = map([mpop.popA, mpop.popB]) do pop
    tsx = Fwd.simplify(ts, pop.nodes)
    rts = Fwd.reverse_relabel(tsx)
    pyts = Fwd.to_tskit(rts)
    xs = collect(pyts.breakpoints())
    ps = pyts.diversity(mode="branch", windows=xs) ./ 2
    xs[2:end], ps
end
P1 = plot(tt[1]...)
hline!([2NA])
P2 = plot(tt[2]...)
hline!([2NB])
plot(P1, P2, size=(500,200))




