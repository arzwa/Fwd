using BenchmarkTools
using Random, Fwd, StatsBase, Distributions


# Benchmark the unlinked homogeneous case.
N = 500
L = 20
s = 0.005
u = s/10
M = Chromosomes(fill(LinearPhysMap(G=1, C=0.0), L))
#M = Chromosomes([Unlinked() for _=1:L])
A = Architecture(fill(BiAllelic(u), L), collect(1:L), M)
Φ = GPMap([HaploidLocus(-s, i) for i=1:L])

rng = Random.seed!(12)
x   = [rand(rng, Bool, L) for _=1:N]
pop = WFPopulation(ploidy=Haploid(), N=N, gpm=Φ, arch=A, x=x)
pop = simulate!(rng, pop, 500)

@code_warntype Fwd.generation!(rng, pop)

# Ideally, we achieve the performance of this:
function unlinked_gen(rng, x, x_, s, u)
    N = length(x)
    L = length(x[1])
    w = map(y->exp(-sum(y)*s), x)
    idx = sample(rng, 1:N, Weights(w), 2N, replace=true)
    for i=1:N
        xi = x[idx[i]]
        xj = x[idx[N+i]]
        for l=1:L
            x_[i][l] = rand(rng) < 0.5 ? xi[l] : xj[l]
        end
    end
    for l=1:L
        n = rand(rng, Poisson(N*u))
        idx = sample(rng, 1:N, n, replace=false)
        for k in 1:n
            x_[idx[k]][l] = !x_[idx[k]][l]
        end
    end
    return x_
end

@code_warntype unlinked_gen(rng, pop.x, pop._x, s, u)

@benchmark unlinked_gen(rng, pop.x, pop._x, s, u)
@benchmark Fwd.generation!(rng, pop)
