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

ngen = 50000
rng = Random.seed!(rand(1:2^32))
pop, qs, ts = let pop=deepcopy(mpop), ts=Fwd.init_ts(mpop, C)
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
    pop, qs, ts
end


P1 = plot(vec(mean(qs, dims=1)[:,:,1]), ylim=(0,1), marker=true, ms=2, color=:black)
vline!([D÷2+0.5], title="\$L=$L, s=$s, m=$m\$", xlabel="deme", ylabel="\$q\$")
x, _, _, d = Fwd.diffdiv(ts, 7, D-7)
P2 = plot(x, d, label="\$T_{7,$(D-7)}\$", color=:gray, alpha=0.7, 
    title="cross-deme coalescence times")
x, _, _, d = Fwd.diffdiv(ts, 1, D)
plot!(x, d, label="\$T_{1,$D}\$", color=:black, legend=:topright)
vline!(xs, label="", color=:salmon, alpha=0.5, lw=2, xlabel="map position")
plot(P1, P2, size=(650,250),margin=5Plots.mm)


# Haploid DMI
L = 2
D = 50
N = 50
m = 0.1
s = 0.1
u = 0.0
C = 0.1
R = LinearMap(C)
xs = [C/2-C/5, C/2+C/5]
pops = map(1:D) do d
    x0 = d > D÷2 ? [[true,false] for _=1:N] : [[false,true] for _=1:N] 
    A = Architecture([BiAllelic(u) for _=1:L], xs)
#    M = GPMap([HaploidTwoLocus(s, s, 0.0, 1, 2)])
    M = GPMap([HaploidTwoLocus(0., 0., -s, 1, 2)])
    pop = WFPopulation(ploidy=Haploid(), N=N, arch=A, 
        gpm=M, recmap=R, x=x0, nodes=collect((d*N+1):(d+1)*N))
end
M = zeros(D,D)
for i=1:D
    i > 1 && (M[i,i-1] = m/2)
    i < D && (M[i,i+1] = m/2)
end
#M[:,1] .= 0.0
#M[:,D] .= 0.0
mpop = MetaPop(pops, M)

ngen = 20000
states = [[0,0],[0,1],[1,0],[1,1]]
rng = Random.seed!(rand(1:2^32))
pop, qs, ts = let pop=deepcopy(mpop), ts=Fwd.init_ts(mpop, C)
    qs = Array{Float64,3}(undef, ngen, D, 4)
    @showprogress for i=1:ngen
        pop = Fwd.generation!(rng, pop, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
        end
        for j=1:D
            pm = proportionmap(pop[j].x)
            qs[i,j,:] .= [haskey(pm, x) ? pm[x] : 0.0 for x in states]
        end
    end
    pop, qs, ts
end

map(1:2250:20000) do t
    p = plot()
    map(1:4) do i
        plot!(vec(qs[t,:,i]), label=join(states[i]), lw=2)
    end
    p
end |> x->plot(x..., legend=:left)

x, _, _, d = Fwd.diffdiv(ts, 7, D-7)
P2 = plot(x, d, label="\$T_{7,$(D-7)}\$", color=:gray, alpha=0.7, 
    title="cross-deme coalescence times")
x, _, _, d = Fwd.diffdiv(ts, 1, D)
plot!(x, d, label="\$T_{1,$D}\$", color=:black, legend=:topright)
vline!(xs, label="", color=:salmon, alpha=0.5, lw=2, xlabel="map position")

