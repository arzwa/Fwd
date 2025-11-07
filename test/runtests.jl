using Fwd
using Test
using Random
using PyCall
using StatsBase
msprime = pyimport("msprime")

function compare_tables(ts1, ts2)
    n1 = ts1.tables.nodes
    n2 = ts2.tables.nodes
    @test all(n1.time .≈ n2.time)  # XXX approx, reverse_relabel prone to float error
    @test all(n1.flags .== n2.flags)
    e1 = ts1.tables.edges
    e2 = ts2.tables.edges
    @test all(e1.left .== e2.left)
    @test all(e1.right .== e2.right)
    @test all(e1.parent .== e2.parent)
    @test all(e1.child .== e2.child)
end

@testset "Round trip `tskit <-> Fwd`" begin
    ts1 = msprime.simulate(10, recombination_rate=1, random_seed=12)
    ts_ = Fwd.from_tskit(ts1)
    ts2 = Fwd.to_tskit(ts_)
    compare_tables(ts1, ts2)
end

@testset "Simplification, test against tskit on msprime sims" begin
    for n in [10, 100, 1000]
        ts0 = msprime.simulate(n, recombination_rate=1, random_seed=1)
        ts1 = Fwd.from_tskit(ts0)
        for N in 2:10
            smple = collect(1:N)
            ts2 = Fwd.simplify(ts1, smple)
            ts3 = Fwd.to_tskit(ts2)
            ts4 = ts0.simplify(smple .- 1)
            compare_tables(ts3, ts4)
            ts5 = Fwd.reverse_relabel(ts1)
            ts6 = Fwd.simplify(ts5, reverse(length(ts5.nodes) .- smple .+ 1))
            ts7 = Fwd.reverse_relabel(ts6)
            ts8 = Fwd.to_tskit(Fwd.reverse_relabel(ts6))
            compare_tables(ts8, ts4)
            #rng = Random.seed!(12)
            #shuffle!(rng, smple)
            #ts9 = Fwd.to_tskit(Fwd.simplify(ts1, smple))
            #compare_tables(ts9, ts4)
        end
    end
end

@testset "Recombination" begin
    # This is a visual test to convince myself of the correctness
    L = 10000
    C = 1.0
    x = fill(true, L)
    y = fill(false, L)
    xs = sort(rand(L)) .* C
    bs = sort(rand(6) .* C)
    z  = similar(x)
    Fwd.recombine!(z, bs, x, y, xs)
    xi = true
    x0 = 1
    for k=1:length(bs)
        i = findfirst(x->x>bs[k], xs)
        @test all(z[x0:i-1] .== xi)
        x0 = i
        xi = !xi
    end
    #@btime Fwd.recombine!(z, bs, x, y, xs)
    #heatmap(xs, 1:3, permutedims([x y z]), size=(500,100))
    #vline!(bs, lw=2, color=:orange, yticks=false, xlim=(0,C), legend=false)
end

@testset "TwoPopOneWay" begin
    NA = 2
    NB = 3
    C = 0.005
    m = 0.2
    A = Architecture(BiAllelic{Float64}[], Float64[], LinearMap(C))
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=A, nodes=collect(1:NA))
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=A, nodes=collect(1:NB) .+ NA)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
    ts = Fwd.init_ts(mpop) 
    rng = Random.seed!(323)
    for i=1:10
        mpop = Fwd.generation!(rng, mpop, ts);
    end
    @test all([ts.nodes[x].pop == 1 for x in mpop.popA.nodes])
    @test all([ts.nodes[x].pop == 2 for x in mpop.popB.nodes])
    mpop, sts = Fwd.simplify!(mpop, ts)
    @test all([sts.nodes[x].pop == 1 for x in mpop.popA.nodes])
    @test all([sts.nodes[x].pop == 2 for x in mpop.popB.nodes])
end

@testset "MetaPop vs. TwoPopOneWay" begin
    NA = 1
    NB = 200
    s = 0.05
    m = s/4
    u = s/200
    AB = Architecture([BiAllelic(  u)], Float64[], Unlinked())
    AA = Architecture([BiAllelic(0.0)], Float64[], Unlinked())
    Φ = GPMap([HaploidLocus(-s, 1)])
    popA = WFPopulation(ploidy=Haploid(), gpm=Φ, N=NA, arch=AA, x=[[true]], nodes=collect(1:NA))
    popB = WFPopulation(ploidy=Haploid(), gpm=Φ, N=NB, arch=AB, x=[[false] for _=1:NB], nodes=collect(1:NB) .+ NA)
    # TwoPop
    twopop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
    rng = Random.seed!(12)
    _,ts1,q1 = simulate!(rng, twopop, init_ts(twopop), 10NB, pop->mean(pop.popB.x)[1])
    # MetaPop
    M = [0.0 m; 0.0 0.0]
    metapop = Fwd.MetaPop([deepcopy(popA), deepcopy(popB)], M)
    rng = Random.seed!(12)
    _,ts2,q2 = simulate!(rng, metapop, init_ts(metapop), 10NB, pop->mean(pop[2].x)[1])
    @test all(q1 .== q2)
    @test all(ts1.edges .== ts2.edges)
end

#@testset "" begin
#    NA = 1000
#    NB = 1000
#    C = 0.5
#    m = 0.001
#    R  = LinearMap(C)
#    nA = collect(1:NA)
#    nB = collect(1:NB) .+ NA
#    A = Architecture(HaploidLocus{Float64}[], Float64[], LinearMap(C))
#    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=A, nodes=collect(1:NA))
#    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=A, nodes=collect(1:NB) .+ NA)
#    mpop = Fwd.TwoPopOneWay(m, popA, popB)
#    ts = Fwd.init_ts(mpop, C) 
#    rng = Random.seed!(152)
#    ngen = 10*(NB+NA)
#    @showprogress for i=1:ngen
#        mpop = Fwd.generation!(rng, mpop, ts);
#        if i % 100 == 0 
#            mpop, ts = Fwd.simplify!(mpop, ts)
#        end
#    end
#    xx, ta, tb, tab = Fwd.diffdiv(ts)
#    Eta = NA
#    Etb = (3NB − 4NB*m + 2NA*NB*m + m^2*NB − m^2*NA*NB)/(1 − 2m + 2NB*m + m^2 − m^2*NB)
#    Etb_ = NB*((3-4m) + m*NA*(2-m))/(1 + 2m*(NB-1))
#    Etb_ = (3 + 2m*NA)/(1 + 2m*NB)
#    Etab = 1/m + NA 
#    @info mean(ta), Eta
#    @info mean(tb), Etb, Etb_*NB
#    @info mean(tab), Etab
#end

