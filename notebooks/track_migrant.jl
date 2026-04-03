
rng = Random.seed!(135)
R = LinearMap(0.1)
AA  = Architecture(BiAllelic[], xs, R)
MA  = GPMap(HaploidLocus[])
NA = 1
NB = 100
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ones(Int, L) for _=1:NA]
xB = [zeros(Int, L) for _=1:NB]
m  = 0.05 
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AA, gpm=MA, x=xB, nodes=nB)
pop = Fwd.TwoPopOneWay(m, popA, popB)

rng = Random.seed!(43)
Fwd.generation!(rng, pop, ts)

ngen = 100
ts = init_ts(pop)
pop, ts, qs_ = Fwd.simulate!(pop, ts, ngen, 
    pop->mean(pop.popB.x), simplify=Inf)

src = filter(k->ts[k].pop == 1, 1:length(ts.nodes))
children = map(src) do n
     children = unique(map(e->e.child, ts.edges[ts.adjlist[n]]))
     filter(x->ts[x].pop == 2, children)
end

nf1 = length.(children)
mean(nf1)
var(nf1)
m*NB

# We seem to have twice the expected variance?
#
# No, not true, this is in fact expected, consider the following
# simulation:

res = map(1:10000) do _
    Z = rand(Poisson(m*NB))  # number of migrants
    Y = mapreduce(_->rand(Poisson(2)), +, 1:Z, init=0)  # number of offspring from migrants (F1s)
    X = rand(Binomial(Y, 0.5))  # number of gene copies from migrant 
end
mean(res), var(res)

# Y has a compound poisson distribution.
