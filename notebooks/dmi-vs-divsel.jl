
using Random, ProgressMeter, Fwd, Distributions, StatsBase

# A haploid DMI is not maintained under mainland-island conditions.
# Architecture:
s = 0.05
m = 0.001
u = s/200
r = m
C = Fwd.distance(r)
AA = Architecture([BiAllelic(0.) for i=1:2], [0.0, C])
AB = Architecture([BiAllelic( u) for i=1:2], [0.0, C])
M = GPMap([HaploidTwoLocus(0.0, 0.0, -s, 1, 2)])
R = LinearMap(C)

states = [[0,0],[0,1],[1,0],[1,1]]
Fwd.eval_component.(Ref(M[1]), states)

# Population model
NA = 1
NB = 500
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [[0,1] for _=1:NA]
xB = [[1,0] for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=M, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=M, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)

states = [[0,0],[0,1],[1,0],[1,1]]
pm = proportionmap(mpop.popB.x)
res = [[haskey(pm, x) ? pm[x] : 0.0 for x in states]]
for i=1:10000
    mpop = Fwd.generation!(Random.default_rng(), mpop)
    pm = proportionmap(mpop.popB.x)
    push!(res, [haskey(pm, x) ? pm[x] : 0.0 for x in states])
end
plot(permutedims(hcat(res...)), 
    label=reshape(join.(states), 1,4), legend=:outertopright)


# Haploid compensatory
s = 0.1
m = 0.001
u = s/500
r = m/100
C = Fwd.distance(r)
AA = Architecture([BiAllelic(0.) for i=1:2], [0.0, C])
AB = Architecture([BiAllelic( u) for i=1:2], [0.0, C])
M = GPMap([Fwd.HaploidTwoLocus(-s, -s, 0.0, 1, 2)])
R = LinearMap(C)
Fwd.eval_component.(Ref(M[1]), states)

NA = 1
NB = 500
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [[1,1] for _=1:NA]
xB = [[0,0] for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=M, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=M, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
states = [[0,0],[0,1],[1,0],[1,1]]
pm = proportionmap(mpop.popB.x)
res = [[haskey(pm, x) ? pm[x] : 0.0 for x in states]]
for i=1:10000
    mpop = Fwd.generation!(Random.default_rng(), mpop)
    pm = proportionmap(mpop.popB.x)
    push!(res, [haskey(pm, x) ? pm[x] : 0.0 for x in states])
end
plot(permutedims(hcat(res...)), 
    label=reshape(join.(states), 1,4), legend=:outertopright)


# Diploid dominant DMI
s = 0.05
m = 0.001
u = s/500
r = m/10
C = Fwd.distance(r)
AA = Architecture([BiAllelic(0.) for i=1:2], [0.0, C])
AB = Architecture([BiAllelic( u) for i=1:2], [0.0, C])
M = GPMap([Fwd.DiploidDominantDMI(-s, 1, 2)])
R = LinearMap(C)

NA = 1
NB = 2000
nA = collect(1:2NA)
nB = collect(1:2NB) .+ 2NA
xA = [[0,1] for _=1:2NA]
xB = [[1,0] for _=1:2NB]
popA = WFPopulation(ploidy=Diploid(), N=NA, arch=AA, gpm=M, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Diploid(), N=NB, arch=AB, gpm=M, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
states = [[0,0],[0,1],[1,0],[1,1]]
pm = proportionmap(mpop.popB.x)
res = [[haskey(pm, x) ? pm[x] : 0.0 for x in states]]
for i=1:10000
    mpop = Fwd.generation!(Random.default_rng(), mpop)
    pm = proportionmap(mpop.popB.x)
    push!(res, [haskey(pm, x) ? pm[x] : 0.0 for x in states])
end
plot(permutedims(hcat(res...)), 
    label=reshape(join.(states), 1,4), legend=:outertopright)


# I would like to compare a model with a bunch of pairwise
# incompatibilities (DMIs) with divergently selected alleles. 
# The relevant comparison is the one where the F1 has the same expected
# fitness.

# DMI model
rng = Random.seed!(67)
L   = 20  # pairs of loci 
s   = 0.3/L
u   = s/200
m   = s/10
C   = 1.0
xs  = collect((C/4L):(C/2L):(C-C/4L))
AA  = Architecture([BiAllelic(0.0) for _=1:2L], xs)
AB  = Architecture([BiAllelic(  u) for _=1:2L], xs)
idx = shuffle(1:2L)
M   = GPMap([HaploidTwoLocus(-s, -s, 0.0, idx[i], idx[i+L]) for i=1:L])
R   = LinearMap(C)
NA  = 1
NB  = 500
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Int, 2L) for _=1:NA]
xB = [zeros(Int, 2L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=M, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=M, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 10^5

#mpop = Fwd.generation!(Random.default_rng(), mpop)


# Simulation
nrep = 1
res  = map(1:nrep) do _
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    pop, ts, qs = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop, C)
        qs=Matrix{Float64}(undef, ngen, 2L)
        @showprogress for i=1:ngen
            pop = Fwd.generation!(rng, pop, ts);
            qs[i,:] .= mean(pop.popB.x)
            if i % 50 == 0 
                pop, ts = Fwd.simplify!(pop, ts)
            end
        end
        pop, ts, qs
    end
    seed, pop, ts, qs
end

plot(res[1][end][1:50000,1:5], ylim=(0,1))

Xs = map(res) do (seed, pop, ts, qs)
    xx, pa, pb, dab = Fwd.diffdiv(ts)
end

xx, yy = Fwd.summarize_wins(first.(Xs), last.(Xs))
p1 = plot(xx, vec(mean(yy, dims=1)), line=:steppre, color=:gray, 
    alpha=0.4, yscale=:log10, ylabel="\$T_{AB}\$")
for l in M.components
    vline!([xs[l.i], xs[l.j]], size=(700,200), color=:lightgray)
end
hline!([1/m + NA])
plot!(xlabel="map position", margin=4Plots.mm, 
    title="$L pair DMIs, $nrep replicates, \$s=$s, m=$m, N=$NA\$")



# With a second continent, the DMI does constitute a barrier I think
s = 0.1
m = 0.005
u = s/1000
r = 10m
C = 0.2 
c = Fwd.distance(r)
xs = [C/2-c/2, C/2+c/2]
AA = Architecture([BiAllelic(0.) for i=1:2], xs)
AB = Architecture([BiAllelic( u) for i=1:2], xs)
M = GPMap([HaploidTwoLocus(-s, 0.0, -s, 1, 2)])
R = LinearMap(C)

states = [[0,0],[0,1],[1,0],[1,1]]
Fwd.eval_component.(Ref(M[1]), states)

# Population model
nrep = 20
res = map(1:nrep) do _
    NA = 1
    NB = 500
    NC = 1
    nA = collect(1:NA)
    nB = collect(1:NB) .+ NA
    nC = collect(1:NC) .+ (NA + NB)
    xA = [[0,1] for _=1:NA]
    xB = [[1,0] for _=1:NB]
    xC = [[1,0] for _=1:NC]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=M, recmap=R, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=M, recmap=R, x=xB, nodes=nB)
    popC = WFPopulation(ploidy=Haploid(), N=NC, arch=AA, gpm=M, recmap=R, x=xC, nodes=nC)
    mpop = MetaPop([popA, popB, popC], [0.0 m 0.0 ; 0.0 0.0 0.0; 0.0 m 0.0])
    ts = init_ts(mpop, C)
    pm = proportionmap(mpop[2].x)
    res = [[haskey(pm, x) ? pm[x] : 0.0 for x in states]]
    @showprogress for i=1:10000
        mpop = Fwd.generation!(Random.default_rng(), mpop, ts)
        if i % 100 == 0
            mpop, ts = Fwd.simplify!(mpop, ts)
        end
        pm = proportionmap(mpop[2].x)
        push!(res, [haskey(pm, x) ? pm[x] : 0.0 for x in states])
    end
    res, ts
end

plot(permutedims(hcat(res[10][1]...)), 
    label=reshape(join.(states), 1,4), legend=:outertopright)

xx = map(last.(res)) do ts
    x1, _, _, dab1 = Fwd.diffdiv(ts, 0, 1)
    x2, _, _, dab2 = Fwd.diffdiv(ts, 2, 1)
    x1, dab1, x2, dab2
end

x1, y1 = Fwd.summarize_wins(getindex.(xx,1), getindex.(xx,2))
plot(x1, vec(mean(y1,dims=1)))
x2, y2 = Fwd.summarize_wins(getindex.(xx,3), getindex.(xx,4))
plot!(x2, vec(mean(y2,dims=1)))
vline!(AA.xs)


