using Fwd, Random, StatsBase, ProgressMeter

function rep(seed, ngen, every)
    NA = 500
    NB = 500
    Ls = 1.0
    L  = 10
    s  = Ls/L
    C  = 1.0
    Nm = 1.0
    m  = Nm/NB
    h  = 0.5
    u  = s*h/200
    xs = collect(C/2L:C/L:C)
    AA = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
    AB = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u  ) for _=1:L], xs)
    R  = LinearMap(C)
    xA = [ ones(Bool, L) for _=1:2NA]
    xB = [zeros(Bool, L) for _=1:2NB]
    nA = collect(1:2NA)
    nB = collect(1:2NB) .+ 2NA
    popA = Fwd.DiploidWFPopulation(N=NA, arch=AA, recmap=R, x=deepcopy(xA), nodes=nA)
    popB = Fwd.DiploidWFPopulation(N=NB, arch=AB, recmap=R, x=deepcopy(xB), nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    ts = Fwd.init_ts(mpop, C) 
    qs  = Vector{Float64}[]
    rng = Random.seed!(seed)
    @showprogress for i=1:ngen
        mpop = Fwd.generation!(rng, mpop, ts)
        push!(qs, mean(mpop.popB.x))
        if i % every == 0
            mpop, ts = Fwd.simplify!(mpop, ts)
        end
    end
    rts = reverse_relabel(ts)
    pts = Fwd.to_tskit(rts)
    mpop, rts, pts, qs
end

res = rep(10, 100*500, 100)

Q = permutedims(hcat(res[end]...))
q̄ = vec(mean(Q, dims=1))

function theights(ts, tmax)
    sts = ts.simplify(ts.samples())
    tms = map(sts.trees()) do tree
        length(tree.roots) > 1 ? tmax : tree.time(tree.root)
    end
    return collect(sts.breakpoints())[2:end-1], tms[1:end-1]
end

x, y = theights(res[3], 100*500)
plot(x,y, size=(900,150), yscale=:log10, color=:salmon, alpha=0.5)
vline!([res[1].popB.arch.xs], color=:black)


