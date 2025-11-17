using Distributed 
addprocs(10)

@everywhere using Pkg
@everywhere Pkg.activate("/home/arzwa/dev/Fwd")
@everywhere using Fwd, Barriers, StatsBase, Distributions, Parameters, ProgressMeter
using Plots; plotsdefault()

@everywhere begin
@with_kw struct Model1{T}
    L :: Int 
    s :: T
    m :: T
    u :: T
    r :: T
    N :: Int
    sb:: T = s
end

function initialp(model)
    @unpack L, s, m, u, r, N, sb = model
    # prediction
    l = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
    A = Barriers.Architecture(l, fill(NaN, L), fill(0.5,L,L))
    M = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
    p = Barriers.Equilibrium(M).Ep[1]
end

function initialize_pop(rng, model, p=initialp(model))
    @unpack L, s, m, u, r, N, sb = model
    # prediction
    # set-up model
    d = Fwd.distance(r)
    lastlocusmap = LinearPhysMap(G=2, C=d)
    otherloci = [LinearPhysMap(G=1, C=0.0) for _=1:L-1]
    cs = Chromosomes([otherloci; lastlocusmap])
    AA = Architecture([BiAllelic(0.0) for _=1:L+1], collect(1:L+1), cs)
    AB = Architecture([[BiAllelic(u) for _=1:L] ; BiAllelic(0.0)], collect(1:L+1), cs)
    ΦA = GPMap([HaploidLocus(0., i) for i=1:L+1])
    ΦB = GPMap([[HaploidLocus(-s, i) for i=1:L]; HaploidLocus(sb, L+1)])
    xA = [[ones( Bool, L); false] for _=1:1]
    xB = [[rand(rng, Bernoulli(1-p), L); false] for _=1:N]
    popA = WFPopulation(ploidy=Haploid(), N=1, gpm=ΦA, arch=AA, x=xA)
    popB = WFPopulation(ploidy=Haploid(), N=N, gpm=ΦB, arch=AB, x=xB)
    pop = TwoPopOneWay(m, popA, popB)
end
end

N  = 500
Ls = 0.5
L  = 20
s  = Ls/L
u  = s/200
m  = s/5

# Unlinked architecture, the last locus has a nonzero-map length as we will
# simulate a linked locally beneficial allele arising there
r = s
d = Fwd.distance(r)
lastlocusmap = LinearPhysMap(G=2, C=d)
otherloci = [LinearPhysMap(G=1, C=0.0) for _=1:L-1]
M = Chromosomes([otherloci; lastlocusmap])
A = Architecture([[BiAllelic(u) for _=1:L] ; BiAllelic(0.0)], collect(1:L+1), M)

function twolocus_cb(pop, i, j; states=[[0,0],[0,1],[1,0],[1,1]])
    pm = proportionmap(map(x->x[[i,j]], pop.x))
    [haskey(pm, x) ? pm[x] : 0.0 for x in states]
end

# Check first if this is simulating linkage correctly
Ds = map(1:20) do _
    Φ = GPMap([HaploidLocus(0., i) for i=1:L+1])
    x = [[ones(Bool, L+1) for _=1:N÷2]; [zeros(Bool, L+1) for _=1:N÷2]]
    pop = WFPopulation(ploidy=Haploid(), N=N, gpm=Φ, arch=A, x=x)
    pop, hs = simulate!(pop, 500, pop->twolocus_cb(pop, L, L+1))
    H = permutedims(hcat(hs...))
    D = H[:,4] .* H[:,1] .- H[:,2] .* H[:,3]
end

plot(mean(Ds), color=:black, ylabel="\$D\$", xlabel="generation")
plot!(x->0.25*exp(-r*x), lw=2)

# Now with selection, check if everything is as expected by comparing
# against mₑ theory predictions.
AA = Architecture([BiAllelic(0.0) for _=1:L+1], collect(1:L+1), M)
AB = Architecture([[BiAllelic(u) for _=1:L] ; BiAllelic(0.0)], collect(1:L+1), M)
ΦA = GPMap([HaploidLocus(0., i) for i=1:L])
ΦB = GPMap([HaploidLocus(-s, i) for i=1:L])
xA = [ones( Bool, L+1) for _=1:1]
xB = [zeros(Bool, L+1) for _=1:N]
popA = WFPopulation(ploidy=Haploid(), N=1, gpm=ΦA, arch=AA, x=xA)
popB = WFPopulation(ploidy=Haploid(), N=N, gpm=ΦB, arch=AB, x=xB)
pop = TwoPopOneWay(m, popA, popB)
cb(pop) = mean(pop.popB.x)

mbs = 0.05:0.1:1.25
qsim = let popA = WFPopulation(ploidy=Haploid(), N=1, gpm=ΦA, arch=AA, x=xA),
           popB = WFPopulation(ploidy=Haploid(), N=N, gpm=ΦB, arch=AB, x=xB)
    qs = map(mbs) do ms
        pop = TwoPopOneWay(ms*s, popA, popB)
        pop, qs = simulate!(pop, 10_000, cb) 
        q̄ = hcat(qs...)[1:end-1,:] |> mean
    end
end

qest = map(range(extrema(mbs)..., 100)) do ms
    l = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
    A = Barriers.Architecture(l, fill(NaN, L), fill(0.5,L,L))
    M = Barriers.MainlandIslandModel(arch=A, m=ms*s, N=N)
    p = Barriers.Equilibrium(M).Ep[1]
    ms, p
end

plot(qest, color=:gray, xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
scatter!(mbs, 1 .- qsim) 



function simulate_fixation!(rng, pop, stopat=1.0)
    L = length(pop.popB.x[1])-1
    k = rand(rng, 1:length(pop.popB.x))
    pop.popB.x[k][L+1] = true
    n = 0
    qs = [mean(pop.popB.x)]
    xs = [twolocus_cb(pop.popB, L, L+1)]
    while true
        n += 1
        pop = generation!(rng, pop)
        q = mean(pop.popB.x)
        push!(qs, q)
        push!(xs, twolocus_cb(pop.popB, L, L+1))
        q[end] == 0.0 && (return (false, qs, xs))
        q[end] >= stopat && (return (true, qs, xs))
    end
end

function sim_nfixations(rng, n, model)
    p = initialp(model)
    pop = initialize_pop(rng, model, p)
    K = 0
    k = 0
    while k < n
        K += 1
        pop = simulate!(rng, pop, model.N) 
        succ, qs = simulate_fixation!(rng, deepcopy(pop))
        succ && (k += 1)
        succ && (@info k, K, k/K)
    end
    return k/K
end

N  = 500
Ls = 0.5
L  = 20
s  = Ls/L
u  = s/200
m  = s/5
r  = s/10
model = Model1(L=L, s=s, m=m, u=u, N=N, r=r, sb=s)
rng = Random.seed!(934)
pop = initialize_pop(rng, Model1(L=L, s=s, m=m, u=u, N=N, r=r, sb=s))
pop, qs = simulate!(pop, N, pop->mean(pop.popB.x)) 
plot(hcat(qs...)'[:,1:L])
# Looks rather equilibrium to me


rng = Random.seed!(932)
pop = initialize_pop(rng, model)
pop = simulate!(pop, N) 

plot()
@showprogress for i=1:50
    _, qs, xs = simulate_fixation!(rng, deepcopy(pop))
    plot!(last.(qs))
end
plot!()



rng = Random.seed!(12)
mm = range(0.05, 1.2, 10)
N  = 500
C  = 0.3
s  = 0.03
LL = [10,20,30]
u  = s/100
r  = 0.03
res = map(LL) do L
    @showprogress map(mm) do ms
        m = ms*s
        K, k = sim_nfixations(rng, 100, 1000, L, s, m, u, N, r, s)
    end
end

plot()
map(zip(Ls, res)) do (ls, Xs)
    p = last.(Xs) ./ first.(Xs)
    plot!(mm, p, marker=true, ms=2)
end
hline!([2sb], color=:gray, ls=:dash)


# Clever (?) MC estimate of fixation probabilities
function initpops(rng, model, n)
    pops = initialize_pops(rng, model, n)
    pops = @showprogress map(pop->simulate!(rng, pop, model.N), pops)
end
function mcfix(rng, model, pops, stopat=0.9999)
    @unpack N, L = model
    n = length(pops)
    # introduce mutation
    for pop in pops
        j = rand(rng, 1:N)
        pop.popB.x[j][L+1] = true
    end
    k = 0
    ps = Float64[]
    while true
        k += 1
        pops′ = map(pop->Fwd.generation!(rng, pop), pops)
        survived = filter(i->sum(last.(pops′[i].popB.x)) > 0, 1:n)
        p = length(survived)/n
        push!(ps, p)
        idx = sample(rng, survived, n)
        pops = [deepcopy(pops′[i]) for i in idx]
        p >= stopat && return ps
    end
end

# check against critical branching process
N = 50
model = Model1(L=L, s=0., m=m, u=u, N=N, r=r, sb=0.)
pops  = initpops(rng, model, 1000)

ps = mcfix(rng, model, deepcopy.(pops), 0.999)
# each gene is expected to have Poisson(1) descendants
# pdf(Poisson(1), 0) = exp(-1)
ps[1], 1-exp(-1)  # 1 - pdf(Poisson(1), 0)
# fixation probability should be 1/N

pp = @showprogress map(1:100) do _
    ps = mcfix(rng, model, deepcopy.(pops), 0.999)
    prod(ps)
end

function mcfix2(rng, model, n; 
        burnin=model.N, stopat=1.0, minX=model.N, p=initialp(model))
    @unpack N, L = model
    K = 0
    pops = @showprogress "first generation" map(1:n) do _
        X = 0
        pop = initialize_pop(rng, model, p)
        pop = simulate!(rng, pop, burnin)
        while true
            K += 1
            j = rand(rng, 1:N)  # individual to mutate 
            pop.popB.x[j][L+1] = true
            pop = Fwd.generation!(rng, pop) 
            X = sum(last.(pop.popB.x))
            X > 0 && break
            # one could reinitialize here? but why not just continue from
            # pop as it is now?
        end
        X, pop
    end
    ps = [n/K]
    while true
        K = 0
        pops = map(1:n) do _
            i = rand(rng, 1:n)  # sample from the set of previous gens
            X, pop = pops[i]    # X is the number of mutant alleles
            K += 1
            # now we need to estimate the probability of survival
            # if X >= minX, we consider survival certain
            # if not, we simulate the next generation
            while true #X < minX  # if X >= minX, we consider it destined to fix
                X >= minX && break
                pop = Fwd.generation!(rng, deepcopy(pop))
                X = sum(last.(pop.popB.x))
                X > 0 && break  # the mutant allele survived
                i = rand(rng, 1:n)
                X, pop = pops[i]
                K += 1
            end
            X, pop
        end
        p = n/K
        push!(ps, p)
        @info prod(ps), p
        p >= stopat && return ps
    end
end

model = Model1(L=0, s=0., m=0., u=0.0, N=200, r=0.499, sb=0.)
rng = Random.seed!(5182)
map(1:10) do _
    res = mcfix2(rng, model, 10000, p=0.5)
end

sim_nfixations(rng, 50, model)

map(1:10) do _
    res = mcfix2(rng, model, 1000, burnin=50, stopat=1.0)
    prod(res)
end

res = map(1:10) do _
    res = mcfix2(rng, model, 1000, burnin=50, stopat=1.0)
    prod(res)
end

sim_nfixations(rng, 1000, 10000, model)

rng = Random.seed!(12)
s  = 0.02
LL = [20,30,40,80]
mm = range(s/10, 3s/4, 10) 
model = Model1(L=10, s=s, m=0.01, u=1e-5, N=500, r=0.05, sb=0.05)
map(LL) do L
    map(mm) do m
        initialp(reconstruct(model, L=L, m=m))
    end 
end

L = 30
s = 0.02
mm = range(s/10, 3s/4, 10)
model = Model1(L=L, s=s, m=0.01, u=1e-5, N=500, r=s, sb=0.05)
res = map(mm) do m
    mcfix2(rng, reconstruct(model, m=m), 10_000, 
        burnin=50, stopat=1.0, minX=100)
end

plot(prod.(res), marker=true, ms=3, color=:black, ylim=(0,0.1))


# --------------------
# This does not work very well...
# I guess resampling (commented in previous fun) does improved things
function mcfix3(rng, model, p0; burnin=model.N, minX=model.N)
    @unpack N, L = model
    pop = initialize_pop(rng, model, p0)
    pop = simulate!(rng, pop, burnin)
    K = 0
    X = 0
    while true
        K += 1
        j = rand(rng, 1:N)  # individual to mutate 
        pop.popB.x[j][L+1] = true
        pop = Fwd.generation!(rng, pop) 
        X = sum(last.(pop.popB.x))
        X > 0 && break
        # one could reinitialize here? but why not just continue from
        # pop as it is now?
    end
    ps = [1/K]
    while true
        K = 0
        while true
            K += 1
            _pop = Fwd.generation!(rng, deepcopy(pop))
            X = sum(last.(_pop.popB.x))
            if X > 0 
                pop = _pop 
                break
            end
        end
        push!(ps, 1/K)
        X >= minX && return ps, exp(sum(log.(ps)))
    end
end

model = Model1(L=L, s=0., m=0., u=1e-5, N=200, r=0.1, sb=0.)
p0 = initialp(model)
res = mcfix3(rng, model, p, burnin=50)

ps = map(1:1000) do _
    _, pfix = mcfix3(rng, model, p0, burnin=50, minX=50)
    pfix
end



# ------------------
# A different take, splitting not be generation, but by mutant frequency
# This seems to work rather well.
@everywhere function mcfix4(rng, model, n; p=initialp(model),
        burnin=model.N, splits=[5,10,20,30,model.N÷2,model.N])
    @unpack L, N = model
    ps = Float64[]
    pops = nothing
    map(enumerate(splits)) do (stage, Xtgt)
        K = 0
        pops = @showprogress "stage $stage" map(1:n) do _
            X, pop = nothing, nothing
            while true
                K += 1
                X, pop = if stage == 1
                    pop = initialize_pop(rng, model, p)
                    pop = simulate!(rng, pop, burnin, show_progress=false)
                    j = rand(rng, 1:N)  # individual to mutate 
                    pop.popB.x[j][L+1] = true
                    1, pop
                else
                    deepcopy(rand(rng, pops))
                end
                while 0 < X < Xtgt
                    pop = Fwd.generation!(rng, pop) 
                    X = sum(last.(pop.popB.x))
                end
                X > 0 && break
            end
            X, pop
        end
        push!(ps, n/K)
        @info (Xtgt, ps[end], prod(ps))
    end
    return ps
end

model = Model1(L=0, s=0., m=0., u=0.0, N=200, r=0.499, sb=0.)
res = mcfix4(rng, model, 10000, p=0.5)

model = Model1(L=20, s=0.05, m=0.02, u=1e-5, N=200, r=0.499, sb=0.05)
res = mcfix4(rng, model, 10000, splits=[5,10,20,40,80,160])


model = Model1(L=50, s=0.02, m=0.005, u=1e-5, N=500, r=0.49, sb=0.05)
@info initialp(model)
rng = Random.seed!(12)
res = mcfix4(rng, model, 1000, splits=[5,10,20,40,80,160])
# 0.0943

model = Model1(L=10, s=0.02, m=0.005, u=1e-5, N=500, r=0.49, sb=0.05)
@info initialp(model)
nrep = 10
seeds = rand(1:2^32, nrep)
pmap(1:nrep) do i
    rng = Random.seed!(seeds[i])
    res = mcfix4(rng, model, 1000, splits=[5,10,20,40,80,160])
end
# 0.085

model = Model1(L=0, s=0.00, m=0.005, u=1e-5, N=500, r=0.49, sb=0.05)
nrep = 10
seeds = rand(1:2^32, nrep)
pmap(1:nrep) do i
    rng = Random.seed!(seeds[i])
    res = mcfix4(rng, model, 1000, splits=[5,10,20,40,80,160], p=0.5)
end
# no divergence -- we expect Pfix ≈ 2(sb - m)


