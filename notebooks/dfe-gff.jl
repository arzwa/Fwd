
using Distributed
@everywhere using Random, Fwd, ProgressMeter, StatsBase, Distributions
using Serialization, Plots; plotsdefault()
using Barriers

rng = Random.seed!(67)
Ls  = 0.5
L   = 25
s̄   = Ls/L
dfe = Exponential(s̄)
u   = 1e-5
m   = 1e-3
C   = 0.5
NA  = 1
NB  = 500
xs   = collect((C/2L):(C/L):(C-C/2L))
ngen = 3*10^5
nrep = 20

res = map(1:nrep) do _
    AA  = Architecture([BiAllelic(0.0) for _=1:L], xs)
    AB  = Architecture([BiAllelic(u) for _=1:L], xs)
    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
    ss  = rand(rng, dfe, L)
    MB  = GPMap([HaploidLocus(-ss[i], i) for i=1:L])
    R   = LinearMap(C)
    nA = collect(1:NA)
    nB = collect(1:NB) .+ NA
    xA = [ ones(Int, L) for _=1:NA]
    xB = [zeros(Int, L) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, 
        gpm=MA, recmap=R, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, 
        gpm=MB, recmap=R, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    seed = rand(rng, 1:2^32)
    hmr = Fwd.hmrecrate(xs)
    (seed, hmr, ss, simulation!(seed, mpop, C, L, ngen)...)
end

xy = map(X->Fwd.diffdiv(X[end-1])[[1,4]], res)

x, Y = Fwd.summarize_wins(first.(xy), last.(xy))

plot(x, vec(mean(Y,dims=1)))

# approximate gff
r̄ = res[1][2]
g = mgf(dfe, -1/r̄)^L
plot(x, vec(mean(Y,dims=1)))
hline!([1/(m*g)])

1 .- map(mean, last.(res))

# Fixed-point iteration using mgf and r̄
function fp(r̄, dfe, m, u, N, L, p=1.0, tol=1e-6)
    s̄ = mean(dfe)
    g = 1.0
    while true
        g = mgf(dfe, -p/(m + r̄))^L
        d = Wright(2N*s̄, N*(m*g + u), N*u, 0.5)
        p_ = 1 - mean(d)
        abs(p - p_) < tol && return p_, g
        p = p_
    end
end

p, g = fp(r̄, dfe, m, u, NB, L)
plot(x, vec(mean(Y,dims=1)))
hline!([1/(m*g)])


