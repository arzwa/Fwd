using Fwd

Ls = [5,10,20,50]
map(Ls) do L
    D = 50
    N = 50
    m = 0.2
    s = 0.01
    u = s/200
    C = 0.1
    R = LinearMap(C)
    xs = collect(C/(2L):C/L:C)
    pops = map(1:D) do d
        sd = d > D÷2 ? s : -s
        x0 = d > D÷2 ? [ones(Bool,L) for _=1:N] : [zeros(Bool,L) for _=1:N] 
        A = Architecture([BiAllelic(u) for _=1:L], xs)
        M = GPMap([HaploidLocus(sd, i) for i=1:L])
        pop = WFPopulation(ploidy=Haploid(), N=N, arch=A, 
            gpm=M, recmap=R, x=x0, nodes=collect((d*N+1):(d+1)*N))
    end
    M = zeros(D,D)
    for i=1:D
        i > 1 && (M[i,i-1] = m/2)
        i < D && (M[i,i+1] = m/2)
    end
    mpop = MetaPop(pops, M)
    ngen = 20000
    rng = Random.seed!(rand(1:2^32))
    pop, qs = let pop=deepcopy(mpop)
        qs = Array{Float64,3}(undef, ngen, D, L)
        @showprogress for i=1:ngen
            pop = Fwd.generation!(rng, pop)
            for j=1:D
                qs[i,j,:] .= mean(pop[j].x)
            end
        end
        pop, qs
    end
    qs
end

plot()
map(zip(Ls,res)) do (L,qs)
    plot!(vec(mean(qs, dims=1)[:,:,1]), label="\$L=$L\$")
end
vline!([D÷2+0.5], legend=:bottomright, label="")


# ts recording
L = 10
D = 20
N = 50
m = 0.2
s = 0.01
u = s/200
C = 0.1
R = LinearMap(C)
xs = collect(C/(2L):C/L:C)
pops = map(1:D) do d
    sd = d > D÷2 ? s : -s
    x0 = d > D÷2 ? [ones(Bool,L) for _=1:N] : [zeros(Bool,L) for _=1:N] 
    A = Architecture([BiAllelic(u) for _=1:L], xs)
    M = GPMap([HaploidLocus(sd, i) for i=1:L])
    pop = WFPopulation(ploidy=Haploid(), N=N, arch=A, 
        gpm=M, recmap=R, x=x0, nodes=collect((d*N+1):(d+1)*N))
end
M = zeros(D,D)
for i=1:D
    i > 1 && (M[i,i-1] = m/2)
    i < D && (M[i,i+1] = m/2)
end
mpop = MetaPop(pops, M)

ngen = 20000
rng = Random.seed!(rand(1:2^32))
pop, qs = let pop=deepcopy(mpop), ts=Fwd.init_ts(mpop, C)
    qs = Array{Float64,3}(undef, ngen, D, L)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
        end
        for j=1:D
            qs[i,j,:] .= mean(pop[j].x)
        end
    end
    pop, qs
end
qs
