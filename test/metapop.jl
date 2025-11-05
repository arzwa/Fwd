using Fwd, Random, WrightDistribution

# special case of mainland-island
N = 500
s = 0.05
u = 0.001

GA = GPMap([HaploidLocus(0., 1)])
GB = GPMap([HaploidLocus(-s, 1)])
AA = Architecture([BiAllelic(0.)], [0.])
AB = Architecture([BiAllelic(u )], [0.])

R = LinearMap(0.)

rng = Random.seed!(123)
pops = [
    WFPopulation(N=1, arch=AA, gpm=GA, recmap=R, ploidy=Haploid(), x=[ ones(Bool,1)]),
    WFPopulation(N=N, arch=AB, gpm=GB, recmap=R, ploidy=Haploid(), x=[zeros(Bool,1) for _=1:N]) 
]
m = s/2
pop = Fwd.MetaPop(pops, [0. m; 0. 0.])

ngen = 10000
qs = Vector{Float64}(undef, ngen)
@showprogress for i=1:ngen
    pop = generation!(rng, pop)
    qs[i] = (sum(pop[2].x) / N)[1]
end

plot(qs, color=:lightgray, alpha=0.5)
#hline!([u/s])
hline!([mean(qs)])
q̄ = mean(Wright(2N*s, N*(m+u), N*u, 0.5))
hline!([q̄])


GA = GPMap([Fwd.DiploidLocus(0., 0., 1)])
GB = GPMap([Fwd.DiploidLocus(-0.5s, -s, 1)])
AA = Architecture([BiAllelic(0.)], [0.])
AB = Architecture([BiAllelic(u )], [0.])
pops = [
    WFPopulation(N=1, arch=AA, gpm=GA, recmap=R, ploidy=Diploid(), x=[ ones(Bool,1) for _=1:2]),
    WFPopulation(N=N, arch=AB, gpm=GB, recmap=R, ploidy=Diploid(), x=[zeros(Bool,1) for _=1:2N]) 
]
pop = MetaPop(pops, [0. m; 0. 0.])
ngen = 10000
qs = Vector{Float64}(undef, ngen)
@showprogress for i=1:ngen
    pop = generation!(rng, pop)
    qs[i] = (sum(pop[2].x) / 2N)[1]
end
