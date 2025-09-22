using Fwd
using Test
using Random
using Parameters
using StatsBase
using WrightDistribution

rng = Random.seed!(12)
N = 1000
L = 1
C = 1e-6
xs= [C/2]
s = 0.05
u = 0.01
A = Architecture([Fwd.HaploidBiLocus(-s, u) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:N]

pop = WFPopulation(N=N, arch=A, recmap=R, ploidy=Haploid(), x=x) 
ts = init_ts(pop, C)
ngen = 10000
qs = Matrix{Float64}(undef, ngen, L)
for i=1:ngen
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


# Poisson locus

rng = Random.seed!(12)
N = 1000
L = 1
C = 1e-6
xs= [C/2]
s = 0.1
u = 0.08
A = Architecture([Fwd.PoissonLocus(-s, u) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Int, L) for _=1:N]

pop = WFPopulation(N=N, arch=A, recmap=R, ploidy=Haploid(), x=x) 
ts = init_ts(pop, C)
ngen = 10000
qs = Matrix{Int}(undef, ngen, N)
for i=1:ngen
    pop = generation!(rng, pop, ts)
    if i % 100 == 0
        pop, ts = Fwd.simplify!(pop, ts)
    end
    qs[i,:] .= first.(pop.x)
end

stephist(vec(qs[ngen÷10:end,:]), bins=0:10, 
    norm=:probability, fill=true, fillalpha=0.2, color=:black)
plot!(0:10, x->pdf(Poisson(u/s),x), line=:steppost)


# BGS
rng = Random.seed!(12)
N = 1000
L = 2
C = 1e-8
xs= [C/2, C/2]
s = 0.08
u = 0.03
un = 0.005
A = Architecture([
    Fwd.PoissonLocus(-s, u), 
    Fwd.HaploidBiLocus(0.0, un)], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [[0, rand() < 0.5 ? 0 : 1] for _=1:N]

pop = WFPopulation(N=N, arch=A, recmap=R, ploidy=Haploid(), x=x) 
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

