using Fwd
using Test
using Random
using Parameters
using StatsBase
using WrightDistribution
using ProgressMeter

rng = Random.seed!(12)
N = 1000
L = 1
C = 1e-6
xs= [C/2]
s = 0.05
u = 0.01
M = GPMap([HaploidLocus(-s, i) for i=1:L])
# To model epistasis in a rather general way, it would be better to revise
# `Architecture` more substantially, think more of it as a
# genotype-phenotyp map. I think one should have each `locus` carry with it
# at which indices in the genotype the relevant alleles are to be found.
# So this becomes more like a quantitative genetics expansion of the
# genotypic value.

A = Architecture([BiAllelic(u) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:N]

pop = WFPopulation(N=N, arch=A, gpm=M, recmap=R, ploidy=Haploid(), x=x) 

ts = init_ts(pop, C)
ngen = 10000
qs = Matrix{Float64}(undef, ngen, L)
@showprogress for i=1:ngen
    pop = generation!(rng, pop, ts)
    if i % 100 == 0
        pop, ts = Fwd.simplify!(pop, ts)
    end
    qs[i,:] .= sum(pop.x) / N
end

plot(qs, color=:lightgray, alpha=0.5)
#hline!([u/s])
hline!([mean(qs)])
q̄ = mean(Wright(2N*s, N*u, N*u, 0.5))
hline!([q̄])


# BGS
rng = Random.seed!(12)
N  = 1000
L  = 2
C  = 1e-8
xs = [C/2, C/2]
s  = 0.08
u  = 0.03
un = 0.005
A  = Architecture([IntAllelic(u), BiAllelic(un)], xs)
M  = GPMap([HaploidLocus(-s, 1), HaploidLocus(0.0, 2)])
R  = LinearMap(C)
ts = Fwd.init_ts(2N, C) 
x  = [[0, rand() < 0.5 ? 0 : 1] for _=1:N]

pop = WFPopulation(N=N, arch=A, gpm=M, recmap=R, ploidy=Haploid(), x=x) 
ts = init_ts(pop, C)
ngen = 50_000
qs = Matrix{Int}(undef, ngen, N)
ps = Vector{Float64}(undef, ngen)
for i=1:ngen
    pop = generation!(rng, pop, ts)
    if i % 100 == 0
        pop, ts = Fwd.simplify!(pop, ts)
    end
    qs[i,:] .= first.(pop.x)
    ps[i] = mean(last.(pop.x))
end

P1 = stephist(vec(qs[ngen÷10:end,:]), bins=0:10, 
    norm=:probability, fill=true, fillalpha=0.2, color=:black)
plot!(0:10, x->pdf(Poisson(u/s),x), line=:steppost)
hline!([exp(-u/s)])
P2 = stephist(ps[ngen÷10:end], norm=true)
Ne = N*exp(-u/s)
d = Wright(0.0, Ne*un, Ne*un, 0.5)
plot!(x->pdf(d, x))
plot(P1, P2, size=(500,200))


# a DMI-like epistatic locus
rng = Random.seed!(12)
N  = 1000
L  = 2
C  = 0.1
xs = [C/3, 2C/3]
s  = 0.08
u  = 0.03
A  = Architecture([BiAllelic(un) for _=1:L], xs)
M  = GPMap([HaploidDMI(-s, 1, 2)])
R  = LinearMap(C)
ts = Fwd.init_ts(2N, C) 
x  = [[true,false] for _=1:N]

function ld(x)
    p1, p2 = mean(x)
    p11 = length(filter(y->y == [1,1], x))/length(x)
    p11 - p1*p2
end

pop = WFPopulation(N=N, arch=A, gpm=M, recmap=R, ploidy=Haploid(), x=x) 
ts = init_ts(pop, C)
ngen = 5000
qs = Matrix{Float64}(undef, ngen, L+1)
for i=1:ngen
    pop = generation!(rng, pop, ts)
    if i % 100 == 0
        pop, ts = Fwd.simplify!(pop, ts)
    end
    qs[i,1:L] .= sum(pop.x) / N
    qs[i,L+1] = ld(pop.x)
end

