using Fwd
using Test
using Random
using PyCall
msprime = pyimport("msprime")

function compare_tables(ts1, ts2)
    n1 = ts1.tables.nodes
    n2 = ts2.tables.nodes
    @test all(n1.time .== n2.time)
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
            @test Fwd.check_edges(ts2, smple)
            ts3 = Fwd.to_tskit(ts2)
            ts4 = ts0.simplify(smple .- 1)
            compare_tables(ts3, ts4)
            ts5 = Fwd.reverse_relabel(ts1)
            ts6 = Fwd.simplify(ts5, length(ts5.nodes) .- smple .+ 1)
            ts7 = Fwd.to_tskit(Fwd.reverse_relabel(ts6))
            compare_tables(ts7, ts4)
        end
    end
end

#@testset "Recombination" begin
#    # This is a visual test to convince myself of the correctness
#    L = 10000
#    C = 1.0
#    x = fill(true, L)
#    y = fill(false, L)
#    xs = sort(rand(L)) .* C
#    bs = sort(rand(6) .* C)
#    z  = similar(x)
#    Fwd.recombine!(z, bs, x, y, xs)
#    #@btime Fwd.recombine!(z, bs, x, y, xs)
#    heatmap(xs, 1:3, permutedims([x y z]), size=(500,100))
#    vline!(bs, lw=2, color=:orange, yticks=false, xlim=(0,C))
#end
#
#@testset "Recombination" begin
#    L = 10000
#    C = 30.0
#    M = LinearMap(C)
#    x = fill(true, L)
#    y = fill(false, L)
#    xs = sort(rand(L)) .* C
#    bs = rand_breakpoints(default_rng(), M)
#    z  = similar(x)
#    Fwd.recombine!(z, bs, x, y, xs)
#    #@btime Fwd.recombine!(z, bs, x, y, xs)
#    heatmap(xs, 1:3, permutedims([x y z]), size=(500,100))
#    vline!(bs, lw=2, color=:orange, yticks=false, xlim=(0,C))
#end
#
#@testset "Haploid fitness" begin
#    L = 1000
#    s = 0.001
#    u = s/100
#    X = rand(Bool, L)
#    x = sort(rand(L))
#    A = Architecture(fill(HaploidBiLocus(-s, u), L), x)
#    @btime fitness(A, X)
#end
#
#@testset "Haploid mutation" begin
#    L = 1000
#    s = 0.001
#    u = s/100
#    X = rand(Bool, L)
#    x = sort(rand(L))
#    A = Architecture(fill(HaploidBiLocus(-s, u), L), x)
#    mutations = Fwd.rand_mutations(rng, A.mut)
#end
#
#@testset "Single locus" begin
#    s = 0.01
#    u = s*0.001
#    m = s*0.1
#    A = Architecture([HaploidBiLocus(-s, u)], [0.0])
#    N = ceil(Int64, 20/s)
#    island  = HaploidWFPopulation([[rand() < m/s] for _=1:N], A, LinearMap(0.))
#    mainland= HaploidFixedPopulation([true]) 
#    metapop = MainlandIsland(island, mainland, m)
#    q = simulate!(metapop, 110N, every=1,
#        callback=P->sum(first.(P.island.x))/P.island.N)[10N:end]
#    d = Wright(2N*s, N*(m + u), N*(u), 0.5)
#    plot(p->pdf(d, p), xlim=(0,1))
#    stephist!(q, bins=0:0.01:1, norm=true)
#    vline!([m/s])
#    vline!([mean(d)])
#    vline!([mean(q)])
#end



