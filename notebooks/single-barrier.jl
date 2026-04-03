using Random, Fwd, ProgressMeter, StatsBase, WrightDistribution
import TreeSequences as TS
using Plots, Serialization; plotsdefault()

NB  = 500
s   = 10/NB
u   = s/200
m   = s/5
C   = 0.1
xs  = [C/2]
L   = 1
R   = LinearMap(C)
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
NA  = 1
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Int, 1) for _=1:NA]
xB = [zeros(Int, 1) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 50NB

# Single example
_, ts, qs = simulate!(rng, deepcopy(mpop), init_ts(mpop), ngen, x->mean(x.popB.x)[1])

stephist(vec(qs), norm=true, bins=0:0.02:1)
plot!(x->pdf(Wright(-2NB*s, NB*u, NB*(u + m), 0.5), 1-x))  

x, ta, tb, tab = TS.diffdiv(ts)
plot(x, tab)

# bunch of replicates
rng = Random.seed!(22)
nrep = 100
res = map(1:nrep) do _
    pop, ts, qs = simulate!(rng, 
        deepcopy(mpop), init_ts(mpop), ngen, x->mean(x.popB.x)[1])
end

tabs = map(res) do (_,ts,_)
    x, ta, tb, tab = Fwd.diffdiv(ts)
end

#serialize("data/tabs-2025-11-24.jls", tabs)
tabs = deserialize("data/tabs-2025-11-24.jls")

q = m/s
p = 1-q
tb_sc(m, q, r) = 1 - 1/r + 1/(q*r) + 1/m

x, tab = Fwd.summarize_wins(first.(tabs), last.(tabs))
plot(x, vec(mean(tab, dims=1)), color=:lightgray)
fr(z) = Fwd.recrate(abs(z-C/2))
plot!(range(extrema(x)..., 200), z->(fr(z) + s)/(m*fr(z)), yscale=:log10, color=:black)
#plot!(range(extrema(x)..., 200), z->(fr(z) + s - m)/(m*fr(z)), yscale=:log10, color=:black, ls=:dash)
plot!(range(extrema(x)..., 200), z->tb_sc(m, q, fr(z)), yscale=:log10, color=:black, ls=:dash)
hline!([ngen])

mₑ(m, r, s) = m*r/(r+s)
tw(NA, NB, m) = NB*(3+2m*NA)/(1+2m*NB)
x, tb = Fwd.summarize_wins(first.(tabs), getindex.(tabs,3))
plot(x, vec(mean(tb, dims=1)), color=:lightgray, legend=:topright, label="")
plot!(range(extrema(x)..., 200), label="\$t_B(m_e)\$",
    x->tw(NA, NB, mₑ(m, fr(x), s*p)))
plot!(range(extrema(x)..., 200), label="hack",
    x->(1-2q)*tw(NA, NB, mₑ(m, fr(x), s*p)) + 2q*(1/mₑ(m, fr(x), s*p) + NA))
plot!(range(extrema(x)..., 200), label="new",
    x->tw2(NA, NB, m, s, fr(x), 1-m/s), yscale=:log10)

function tw2(NA, NB, m, s, r, p)
    me = mₑ(m, r, s*p)
    _tab = 1/me + NA
    _tb = tw(NA, NB, me)
    x = r + s*p
    2m/x^2 + 2m/x * _tab + (x-2m)/x * _tb
end

plot(range(extrema(x)..., 200), x->m/(s*(1-q) + fr(x)))

plot!(range(extrema(x)..., 200), z->tw_sc2(m, q, fr(z), NB), color=:black, ls=:dash)
plot!(range(extrema(x)..., 200), z->tw_sc3(m, q, fr(z), NB), color=:red, ls=:dash)
plot!(ylim=(0,ngen))

function gettabat(tabs, x)
    map(tabs) do (xs,tab)
        j = findfirst(i->xs[i] < x <= xs[i+1], 1:length(xs)-1)
        tab[j]
    end
end

y = 0.02
rr = Fwd.recrate(xs[1] - y)
tx = gettabat(tabs, y)
mean(tx)


(-2*N^2*m^2*q^3*r + 4*N^2*m^2*q^2*r^2 - 2*N^2*m^2*q^2*r - 4*N^2*m^2*q*r^2 + 10*N^2*m^2*q*r - 6*N^2*m^2*r + 4*N^2*m*q^3*r^3 - 6*N^2*m*q^3*r^2 - 4*N^2*m*q^2*r^3 + 18*N^2*m*q^2*r^2 - 12*N^2*m*q*r^2 + 6*N^2*q^3*r^3 - 6*N^2*q^2*r^3 + 2*N*m^2*q^2*r - 4*N*m^2*q*r + 4*N*m^2*q - 4*N*m^2 + 2*N*m*q^3*r^2 - 2*N*m*q^2*r^2 + 3*N*m*q^2*r - 2*N*m*q*r^2 - 3*N*m*q*r - 3*N*m*r + 3*N*q^3*r^2 - 3*N*q^2*r^2 - 3*N*q*r^2 + 2*m*q^2*r - 2*m*q*r + 2*m*q - 2*m + 2*q^3*r - 2*q*r)/(r*(4*N^2*m^2*q^2*r - 4*N^2*m^2*q*r + 4*N^2*m*q^3*r^2 - 4*N^2*m*q^2*r^2 - 2*N*m^2 + 4*N*m*q^2*r - 6*N*m*q*r + 2*N*q^3*r^2 - 2*N*q^2*r^2 - m - q*r))


# With MetaPop
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
mpop = MetaPop([popA, popB], [0.0 m; 0.0 0.0])
ngen = 10^5

rng = Random.seed!(22)
pop2, ts2, qs2 = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop)
    qs = Matrix{Float64}(undef, ngen, 1)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        qs[i,:] .= mean(pop[2].x)
        if i % 50 == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    pop, ts, qs
end

stephist(vec(qs2), norm=true, bins=0:0.02:1)
plot!(x->pdf(Wright(-2NB*s, NB*u, NB*(u + m), 0.5), 1-x))  

x, ta, tb, tab = Fwd.diffdiv(ts)
plot(x, tab)
x, ta, tb, tab = Fwd.diffdiv(ts2)
plot!(x, tab)

# Check the simplified ts  
@info length(filter(x->x.time == 0, ts.nodes))
@info length(filter(x->x.time == ngen, ts.nodes))

# Additional generations without ts simplification
pop, ts = let ngen=1000, pop=deepcopy(pop), ts=deepcopy(ts)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
    end
    pop, ts
end

ns = filter(i->ts[i].time > ngen, 1:length(ts.nodes))
f1s = zeros(1001)
for n in ns
    parpop = ts[n].pop
    childnodes = unique([ts.edges[e].child for e in ts.children[n]])
    childpops = Int[ts[c].pop for c in childnodes] 
    f1s[ts[n].time+1-ngen] += sum(parpop .!= childpops)
end
plot(f1s)
hline!([m*NB])
hline!([mean(f1s)])
    
# introduce a migrant
let dpop = reconstruct(deepcopy(pop), m=0.0), dts = deepcopy(ts)
    k = rand(rng, 1:dpop.popB.N)
    nk = Fwd.migrate!(dpop.popA, dpop.popB, 1, k)
    W = Fwd.eval_fitness(dpop.popB)
    W[k]
end

function trace_descendants(ts, node, C, tmax)
    segments = Tuple{Int,Int,Float64,Float64}[]
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
                push!(segments, (time(ts[edge.child]) - t, 
                    edge.child, z0, z1))
                _recurse(edge.child, z0, z1)
            end
        end
    end
    _recurse(node, 0.0, C)
    segments
end

nrep = 1000
xxs = 0.001:0.001:C-0.001
res = map(1:nrep) do _
    let dpop = reconstruct(deepcopy(pop), m=0.0), dts = deepcopy(ts)
        k = rand(rng, 1:dpop.popB.N)
        nk = Fwd.migrate!(dpop.popA, dpop.popB, 1, k)
        W = Fwd.eval_fitness(dpop.popB)
        for i=1:10
            dpop = Fwd.generation!(rng, dpop, dts)
        end
        ds = trace_descendants(dts, nk, C, 10)
        map(x->length(filter(d->d[2] <= x <= d[3], ds)), xxs)
    end
end
plot(xxs, mean(res))

nrep = 1000
xxs = 0.001:0.001:C-0.001
res = @showprogress map(1:nrep) do _
    let dpop=reconstruct(deepcopy(pop), m=0.0), dts=deepcopy(ts), T=1
        k = rand(rng, 1:dpop.popB.N)
        nk = Fwd.migrate!(dpop.popA, dpop.popB, 1, k)
        residents = filter(x->x!=nk, dpop.popB.nodes) 
        W = Fwd.eval_fitness(dpop.popB)
        for i=1:T
            dpop = Fwd.generation!(rng, dpop, dts)
        end
        rs = map(nx->trace_descendants(dts, nx, C, T), residents)
        ds = trace_descendants(dts, nk, C, T)
        wmig = map(x->length(filter(d->d[2] <= x < d[3], ds)), xxs)
        wres = map(rs) do resdes
            map(x->length(filter(d->d[2] <= x < d[3], resdes)), xxs)
        end
        wmig ./ mean(wres)
    end
end

function backcross_fitnesses(rng, pop, ts, T, x)
    # introduce migrant
    k = rand(rng, 1:pop.popB.N)
    j = rand(rng, 1:pop.popA.N)
    migrant = Fwd.migrate!(pop.popA, pop.popB, j, k)
    residents = filter(x->x!=migrant, pop.popB.nodes) 
    # do T generations
    map(1:T) do t
        pop = Fwd.generation!(rng, pop, ts)
    end
    mdesc = trace_descendants(ts, migrant, C, T)
    rdesc = map(r->trace_descendants(ts, r, C, T), residents)
    map(1:T-1) do t
        # descendants in generation t that inherited from n at x
        dtm = filter(y->y[1] == t && y[3] <= x < y[4], mdesc)
        # those guys' descendants in the next generation
        mdescxt = map(d->length(trace_descendants(ts, d[2], C, 1)), dtm)
        rdescxt = map(rdesc) do r
            dtr = filter(y->y[1] == t && y[3] <= x < y[4], r)
            rdescxt = map(d->length(trace_descendants(ts, d[2], C, 1)), dtr)
            length(dtr) == 0 ? 0.0 : mean(rdescxt)
        end 
        length(mdescxt) == 0 ? 0.0 : mean(mdescxt) / mean(rdescxt)
    end
end
    
dd = backcross_fitnesses(rng,
    reconstruct(deepcopy(pop), m=0.0), 
    deepcopy(ts), 
    5, xs[1])

wm = map(1:1000) do _
    backcross_fitnesses(rng,
        reconstruct(deepcopy(pop), m=0.0), 
        deepcopy(ts), 
        1, xxs)
#    length(md) / mean(length.(rd))
end

mean(wm)

nrep = 1000
ww = map(1:nrep) do _
    dd = backcross_fitnesses(rng,
        reconstruct(deepcopy(pop), m=0.0), 
        deepcopy(ts), 
        1, xxs)
    length(dd)
end |> mean

# dd[x][t] contains for each t-gen descendants that inherited from the
# migrant at x the number of ofsspring in the next generation

nrep = 1000
xxs = xs
ww = map(1:nrep) do _
    dd = backcross_fitnesses(rng,
        reconstruct(deepcopy(pop), m=0.0), 
        deepcopy(ts), 
        1, xxs)
    length.(last.(last.(first.(dd))))
end;
mean(ww)

w2 = map(enumerate(xxs)) do (i,x)
    zs = getindex.(ww, i)
    map(1:5) do t
        mean(vcat(getindex.(zs, t)...))
    end
end

Wm = hcat(mean(ww)...)
plot(xxs, Wm')



# ----------------------------------------------------------------------
# Haplotype blocks
NB  = 500
s   = 10/NB
u   = s/200
m   = s/5
C   = 0.1
xs  = [C/2]
L   = 1
R   = LinearMap(C)
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
NA  = 1
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Int, 1) for _=1:NA]
xB = [zeros(Int, 1) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 50NB

rng = Random.seed!(13)
pop, ts, qs = simulate!(rng, 
    deepcopy(mpop), init_ts(mpop), ngen, x->mean(x.popB.x)[1])

stephist(qs)
       
ts = Fwd._add_grand_ancestor(ts)
x, ta, tb, tab = Fwd.diffdiv(ts)
plot(x, tab, yscale=:log10)
plot!(x, tb)

pts = to_tskit(ts)
i = pts.samples(population=0)[1]
j, k = sample(pts.samples(population=1), 2, replace=false)

sts = pts.simplify([i,j,k])
sx = collect(sts.breakpoints())
tres = [x.copy() for x in sts.trees()][1:end-1]

using NewickTree, SmoothTrex
nts = map(x->readnw(x.newick()), tres)

# three possible trees, defined by 
# (1,2) = (B1,B2)
# (1,3) = (B1,A)
# (2,3) = (B2,A)

tops = map(topologize, nts)

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

plot(x, tab, size=(700,200), yscale=:log10)
plot!(x, tb)
for (j,cl) in enumerate([(1,2), (1,3), (2,3)])
    for i=1:length(nts)
        tree = nts[i]
        x0 = sx[i]
        x1 = sx[i+1]
        node = getlca(tree, string.(cl)...)
        top1 = length(getleaves(node)) == 2 
        if top1 
            t1 = distance(node[1])
            t2 = distance(node)
            plot!(rectangle(x1-x0, t2, x0, t1), color=j, line=false, alpha=0.5)
        end
    end
end
plot!()

