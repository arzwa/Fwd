
NA = 20
NB = 10
L = 2
C = 0.1
s = 0.2
m = 0.1
h = 0.5
u = s*h/200
xs = collect(C/2L:C/L:C)
AA = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
AB = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u  ) for _=1:L], xs)
R  = LinearMap(C)
    
ts = Fwd.init_ts(2NA + 2NB, C) 
xA = [ones(Bool, L) for _=1:2NA]
xB = [zeros(Bool, L) for _=1:2NB]
nA = collect(1:2NA)
nB = collect(1:2NB) .+ 2NA
popA = Fwd.DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=xA, nodes=nA)
popB = Fwd.DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)

rng = Random.seed!(322)
mpop = Fwd.generation!(rng, mpop, ts)
Fwd.check_edges(ts, Fwd.active_nodes(mpop))

_popA, tsa = Fwd.simplify!(mpop.popA, ts)
Fwd.check_edges(tsa, _popA.nodes)

mpop, ts = Fwd.simplify!(mpop, ts)
Fwd.check_edges(ts, Fwd.active_nodes(mpop))
 
nodes = Fwd.active_nodes(mpop)
es = [(n, filter(x->x.child == n, ts.edges)) for n in nodes]
ls = map(ex->sum([x.rght-x.left for x in ex[2]]), es)

filter(e->e.parent ∈ nA && e.child ∈ mpop.popB.nodes, ts.edges)

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




