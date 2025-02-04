using Fwd, Plots, Test, BenchmarkTools, WrightDistribution
default(grid=false, fontfamily="Computer modern", framestyle=:box, legend=false)

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
    #@btime Fwd.recombine!(z, bs, x, y, xs)
    heatmap(xs, 1:3, permutedims([x y z]), size=(500,100))
    vline!(bs, lw=2, color=:orange, yticks=false, xlim=(0,C))
end

@testset "Recombination" begin
    L = 10000
    C = 30.0
    M = LinearMap(C)
    x = fill(true, L)
    y = fill(false, L)
    xs = sort(rand(L)) .* C
    bs = rand_breakpoints(default_rng(), M)
    z  = similar(x)
    Fwd.recombine!(z, bs, x, y, xs)
    #@btime Fwd.recombine!(z, bs, x, y, xs)
    heatmap(xs, 1:3, permutedims([x y z]), size=(500,100))
    vline!(bs, lw=2, color=:orange, yticks=false, xlim=(0,C))
end

@testset "Haploid fitness" begin
    L = 1000
    s = 0.001
    u = s/100
    X = rand(Bool, L)
    x = sort(rand(L))
    A = Architecture(fill(HaploidBiLocus(-s, u), L), x)
    @btime fitness(A, X)
end

@testset "Haploid mutation" begin
    L = 1000
    s = 0.001
    u = s/100
    X = rand(Bool, L)
    x = sort(rand(L))
    A = Architecture(fill(HaploidBiLocus(-s, u), L), x)
    mutations = Fwd.rand_mutations(rng, A.mut)
end

@testset "Single locus" begin
    s = 0.01
    u = s*0.001
    m = s*0.1
    A = Architecture([HaploidBiLocus(-s, u)], [0.0])
    N = ceil(Int64, 20/s)
    island  = HaploidWFPopulation([[rand() < m/s] for _=1:N], A, LinearMap(0.))
    mainland= HaploidFixedPopulation([true]) 
    metapop = MainlandIsland(island, mainland, m)
    q = simulate!(metapop, 110N, every=1,
        callback=P->sum(first.(P.island.x))/P.island.N)[10N:end]
    d = Wright(2N*s, N*(m + u), N*(u), 0.5)
    plot(p->pdf(d, p), xlim=(0,1))
    stephist!(q, bins=0:0.01:1, norm=true)
    vline!([m/s])
    vline!([mean(d)])
    vline!([mean(q)])
end



