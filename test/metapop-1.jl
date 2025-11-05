
# Comparison of MetaPop against TwoPopOneWay
using Random, Fwd, ProgressMeter, StatsBase, WrightDistribution
using Plots; plotsdefault()

rng = Random.seed!(67)
s   = 0.05
u   = s/200
m   = s/10
C   = 0.1
xs  = [C/2]
L   = 1
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs)
AB  = Architecture([BiAllelic(  u) for _=1:L], xs)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus( -s, i) for i=1:L])
R   = LinearMap(C)
NA  = 1
NB  = 500
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Int, 1) for _=1:NA]
xB = [zeros(Int, 1) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 10^5

ts = Fwd.init_ts(mpop,C)
rng = Random.seed!(22)
@btime Fwd.generation!(rng, mpop);
@code_warntype Fwd.generation!(rng, mpop);
#    32.887 μs (2139 allocations: 90.68 KiB)


rng = Random.seed!(22)
pop, ts, qs = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop, C)
    qs = Matrix{Float64}(undef, ngen, 1)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        qs[i,:] .= mean(pop.popB.x)
        if i % 50 == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    pop, ts, qs
end

stephist(vec(qs), norm=true, bins=0:0.02:1)
plot!(x->pdf(Wright(-2NB*s, NB*u, NB*(u + m), 0.5), 1-x))  

x, ta, tb, tab = Fwd.diffdiv(ts)
plot(x, tab)

# With MetaPop
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, recmap=R, x=xB, nodes=nB)
mpop = MetaPop([popA, popB], [0.0 m; 0.0 0.0])
ngen = 10^5

ts = Fwd.init_ts(mpop,C)
rng = Random.seed!(22)
@btime Fwd.generation!(rng, mpop);
@code_warntype Fwd.generation!(rng, mpop);
#  32.486 μs (2138 allocations: 89.09 KiB)

rng = Random.seed!(22)
pop2, ts2, qs2 = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop, C)
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

stephist(vec(qs), norm=true, bins=0:0.02:1)
stephist!(vec(qs2), norm=true, bins=0:0.02:1, ls=:dash)
plot!(x->pdf(Wright(-2NB*s, NB*u, NB*(u + m), 0.5), 1-x))  

x2, _, _, tab2 = Fwd.diffdiv(ts2)
plot(x, tab)
plot!(x2, tab2, ls=:dash)
