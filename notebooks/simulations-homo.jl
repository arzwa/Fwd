
# Simulation experiment to assess predictions for homogeneous
# architectures. 
using Serialization
using DataFrames

@everywhere function equallyspaced(r)
    d  = Fwd.distance(r)
    xs = cumsum(fill(d, L))
    C  = last(xs) + d
    R  = LinearMap(C)
    xs, R
end

@everywhere function getmodel(L, m, s, r, N, u)
    xs, R = equallyspaced(r)
    AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
    AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
    MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
    MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
    NA  = 1
    nA  = collect(1:NA)
    nB  = collect(1:N ) .+ NA
    xA  = [ ones(Bool, L) for _=1:NA]
    xB  = [zeros(Bool, L) for _=1:N ]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=N , arch=AB, gpm=MB, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
end

L   = 20 
s   = 0.1/L
Ns  = 5.
N   = ceil(Int, Ns/s)
u   = s/200
rss = [0.2, 1.0, 5.0]
res = map(rss) do rs
    xs, R = equallyspaced(rs*s)
    ms = 1.0 
    # relevant ms range
    while true
        pex = Equilibrium(BPModel(m=ms*s, xs=xs, s=fill(s, L), u=u, Ne=N), α=0.1).Ep[L÷2]
        pex < 0.05 && break
        ms += 0.3
    end
    ms_max = ms + 0.1
    @info rss, ms_max
    mss = range(0.05, ms_max, 15)
    reps = map(1:5) do _
        res = pmap(mss) do ms
            model = getmodel(L, ms*s, s, r, N, u)
            ngen = 2000N
            evry = ceil(Int,N/10)
            @info ms, ngen, evry
            pop, qs = simulate!(model, ngen, x->mean(x.popB.x), every=evry)
            qs
        end
    end
    mss, reps
end

plot()
map(res) do (mss, reps)
    P = map(reps) do qs
        mean(pcat(qs...), dims=2)
    end |> mean |> x->pcat(x...)
    scatter!(mss, 1 .- P[:,10])
end
plot!()

mss = range(0.05, 2, 15)

reps = map(1:5) do _
    res = pmap(mss) do ms
        model = getmodel(L, ms*s, s, r, N, u)
        ngen = 1000N
        evry = ceil(Int,N/10)
        @info ms, ngen, evry
        pop, qs = simulate!(model, ngen, x->mean(x.popB.x), every=evry)
    end
end

P = map(reps) do res
    map(last.(res)) do qs
        mean(pcat(qs...)[end-1000:end,:], dims=1) |> vec
    end  |> x->1 .- hcat(x...)' 
end 

PM = mean(P)

@unpack L, s, r, N, u, ms, P = data[1]
mss2 = range(extrema(mss)..., 100)
PBP = map(mss2) do ms
    Equilibrium(BPModel(m=ms*s, xs=xs, s=fill(s, L), u=u, Ne=N), α=0.1).Ep
end |> x->pcat(x...)

PRV = map(mss2) do ms
    loci = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
    A = Barriers.Architecture(loci, xs, Fwd.rec_matrix(xs))
    M = Equilibrium(Barriers.MainlandIslandModel(arch=A, 
        m=ms*s, N=N, mode=1))
    M.Ep
end |> x->pcat(x...)

# add to results
data = deserialize("data/homo-sims.jls")
serialize("data/homo-sims.jls", data)

map(data) do res
    @unpack BP, PRV, P = res
    PM = mean(res.P)
    plot(BP[1], BP[2][:,10], color=1)
    plot!(PRV[1], PRV[2][:,10], color=2)
    scatter!(mss, PM[:,10], color=:black, ms=2)
    plot!()
end |> x->plot(x...)

