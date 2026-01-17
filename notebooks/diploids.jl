
using Fwd
using Test
using Random
using Parameters
using StatsBase
using ProgressMeter
using Plots
using PyCall
using QuadGK
msprime = pyimport("msprime")

function LDxy(x, y)
    @assert length(x) == length(y)
    n = length(x)  # number of individuals
    ld = 0.0
    k = 0
    for i=1:n-1
        for j=i+1:n
            ld += x[i]*x[j]*y[i]*y[j]
            k += 1
        end
    end
    ld*2/(n*(n-1))
end

function LDxy(G::Matrix, x)
    n = length(x)
    ds = Tuple{Float64,Float64}[]
    for i=1:n-1
        for j=i+1:n
            ld = LDxy(G[i,:], G[j,:])
            d = x[j] - x[i]
            push!(ds, (d, ld))
        end
    end
    return ds
end

function getbins(xy, bins)
    n = length(bins)
    zs = zeros(n-1)
    ks = zeros(Int, n-1)
    a, b = extrema(bins)
    for (x,y) in xy
        i = ceil(Int, (x / (b-a)) * (n-1))
        zs[i] += y
        ks[i] += 1
    end
    bins, zs ./ ks
end

pheal(λ, t) = 0.5 + (exp(-2λ*t)-1)/(4λ*t)
psurvival(λ, x) = quadgk(t->exp(-2*t*x*(1-pheal(λ, t)))*λ*exp(-λ*t), 0, Inf)[1]

function predictbins(bins, λ)
    n = length(bins)
    map(2:n) do i
        quadgk(x->psurvival(λ,x), bins[i-1], bins[i])[1]/(bins[i]-bins[i-1])
    end
end

rng = Random.seed!(8)
nrep = 10
N = 500
C = 0.1
A = Architecture(recmap=LinearMap(C))
res = map(1:nrep) do it
    pop = WFPopulation(N=N, arch=A, ploidy=Diploid()) 
    pop, ts = simulate!(pop, init_ts(pop), 10N)
    pts = to_tskit(ts)
    #plot(Fwd.theights(pts))
    # take a sample
    n = 100
    idx = sample(rng, 1:N, n, replace=false)
    ind = [x for x in pts.samples()]
    smpl = [ind[idx]; ind[idx .+ N]]
    sts = pts.simplify(smpl)
    # simulate mutations
    mts = msprime.sim_mutations(
        sts, rate=1, 
        random_seed=rand(rng, 1:2^32), 
        model=msprime.BinaryMutationModel(), 
        discrete_genome=false)
    x = map(v->v.position, mts.variants())
    H = mts.genotype_matrix()
    G = H[:,1:n] .+ H[:,n+1:end]
    p = mean(H, dims=2)
    GG = permutedims(
        mapreduce(i->(G[i,:] .- 2p[i]) ./ √(2p[i]*(1-p[i])), hcat, 1:length(p)))
    maf = filter(i->0.25 < p[i] < 0.75, 1:length(p))
    GG = GG[maf,:]
    x  = x[maf,:]
    lds = LDxy(GG, x)
#    res = map(lds) do (d, ld)
#        psurvival(1/2N, d), ld
#    end
end

#P0 = scatter(sample(res, 1000, replace=false), ms=2, color=:lightgray, title="\$N_e = $N\$, 0.5M chromosome, \$p > 0.25\$")
#plot!(x->x, lw=2, color=:black, ylabel="\$\\mathrm{LD}_{x,y}\$", xlabel="\$S(u)\$")
bins = 0:0.005:C
zs = map(res) do lds
    _, zs = getbins(lds, bins)
    zs
end
bm = [(bins[i]+bins[i-1])/2 for i=2:length(bins)]

ss = predictbins(bs, 1/N)  # XXX coal rate should be 1/2N ? 
plot(bm, ss, yscale=:log10, marker=false, xticks=[0.02,0.04,0.06,0.08,0.1])
P1 = scatter!(bm, mean(zs), yscale=:log10, color=:black, ms=3,
    xlabel="distance (M)", label="\$\\mathrm{LD}_{x,y}\$")

plot!(bs[2:end], ss, lw=1, color=:gray, label="\$\\overline{S(u)}\$", legend=:topright)

plot(P0,P1,size=(500,200),margin=2Plots.mm)

