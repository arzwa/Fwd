
rng = Random.seed!(29)

nrep = 100
res = map(1:nrep) do k
    N = 500
    C = 0.1
    xs= [C/2]
    s = 0.02
    u = 0.001
    M = GPMap([HaploidLocus(s, 1)])
    R = LinearMap(C)
    A = Architecture([BiAllelic(0.)], xs, R)
    x = [[false] for _=1:N]
    pop = WFPopulation(N=N, arch=A, gpm=M, ploidy=Haploid(), x=x, nodes=collect(1:N)) 
    ngen = 2000
    pop, ts, qs = simulate!(pop, init_ts(pop), ngen, x->mean(x.x), every=1)
    _pop = deepcopy(pop)
    _ts = deepcopy(ts)
    pop.x[rand(1:N)][1] = true
    _qs = [[1/N]]
    while true
        while 0 < mean(pop.x)[1] < 1
            pop = generation!(pop, ts)
            push!(_qs, mean(pop.x))
        end
        mean(pop.x)[1] == 1 && break
        pop = deepcopy(_pop)
        ts = deepcopy(_ts)
        _qs = [[1/N]]
        pop.x[rand(1:N)][1] = true
    end
    ts = TS.simplify(ts, pop.nodes)
    pop, ts, [qs; _qs]
end

t = map(res) do (pop, ts, qs)
    x, y = TS.div(TS._add_grand_ancestor(ts))
end

x, y = Fwd.summarize_wins(t)

plot(x, mean(y, dims=1)' , color=:black, size=(350,220),
    ylabel="\$\\overline{T}_2\$", xlabel="map position (M)")
hline!([N], lw=2, color=:firebrick)
annotate!(0.105, N, text("\$N=500\$", :left, 10, :firebrick),  
    right_margin=15Plots.mm)

k = 25
pts = TS.to_tskit(res[k][2]).simplify(1:10)
pp = pts.diversity(pts.samples(), mode="branch", windows=collect(pts.breakpoints()))
bps = collect(pts.breakpoints())
ts = [t.copy() for t in pts.trees()]
nt = length(ts)
#P0 = plot(t[k], color=:black)
P0 = plot(bps, [pp[1] ; pp], line=:steppost, color=:black, ylabel="\$T\$")
ymn, ymx = ylims(P0)
k1, k2, k3 = 9, nt÷2+6, nt-18
for j = [k1, k2, k3]
    plot!(P0, rectangle(bps[j], bps[j+1], 
        0, ymx), fill=true, lw=0,
        color=:gray, alpha=0.3, ylim=(ymn,ymx))
end
P1 = plot(readnw(ts[k1].as_newick()), orientation=4, ytick=true,
    yshowaxis=true, framestyle=:default, ylabel="\$t\$ before present")
bps
h1 = (-)(ylims(P1)...)
P2 = plot(readnw(ts[k2].as_newick()), orientation=4)
h2 = (-)(ylims(P2)...)
P3 = plot(readnw(ts[k3].as_newick()), orientation=4)
h3 = (-)(ylims(P3)...)
H = abs(min(h1,h2,h3))*1.05
ps = map([P1, P2, P3]) do p
    y0, y1 = ylims(p)
    plot(p, ylims=(y0, y1 + H + y0), yticks=([y0, y0 + 1000, y0 + 2000], [0, 1000, 2000 ]))
end 
Pb = plot(ps..., layout=(1,3))
plot(P0, Pb, layout=(2,1))





