using Random, Fwd, ProgressMeter

rng = Random.seed!(67)
L   = 10   
s   = 0.25/L
u   = s/200
m   = s/10
C   = 1.0
xs  = collect((C/2L):(C/L):(C-C/2L))
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs)
AB  = Architecture([BiAllelic(u) for _=1:L], xs)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
R   = LinearMap(C)
NA  = 500
NB  = 500
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Int, L) for _=1:NA]
xB = [zeros(Int, L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, recmap=R, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, recmap=R, x=xB, nodes=nB)
mpop = Fwd.TwoPopOneWay(m, popA, popB)
ngen = 10^5

# Simulation
nrep = 1
res  = map(1:nrep) do _
    seed = rand(1:2^32)
    rng = Random.seed!(seed)
    pop, ts, qs = let pop=deepcopy(mpop), ts=Fwd.init_ts(pop, C)
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
    seed, pop, ts, qs
end

pop = res[1][2]
ts = res[1][3]

na = nb = 10
smpl = sample(rng, pop.popA.nodes, na, replace=false)
smpl = [smpl; sample(rng, pop.popB.nodes, na, replace=false)]

tss = Fwd.simplify(ts, smpl)


