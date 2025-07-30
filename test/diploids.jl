using Fwd
using Test
using Random
using Parameters
using StatsBase

N = 150
L = 0
C = 1e-12
xs= Float64[]
A = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:2N]

rng = Random.seed!(12)
pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
ws  = Fwd.eval_fitness(pop)
for _=1:10
    idx = sample(rng, 1:N, Weights(ws), 2N)
    pop = Fwd.generation!(rng, pop, idx, ts)
end

rng = Random.seed!(57)
rts = map(1:100) do i
    ts  = Fwd.init_ts(2N, C) 
    x   = [zeros(Bool, L) for _=1:2N]
    pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    @info i
    for i=1:5000
        idx = sample(rng, 1:N, 2N)
        pop = Fwd.generation!(rng, pop, idx, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    Fwd.reverse_relabel(ts)
end
mean(x->time(x.nodes[end]), rts), 4N*(1-1/2N)


# -----------------------------------------
N = 3
L = 0
C = 0.1
xs= Float64[]
A = Architecture([Fwd.DiploidBiLocus(0.0, 0.0, 0.0) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:2N]
pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
rng = Random.seed!(19)
for _=1:3
    idx = sample(rng, 1:N, 2N)
    pop = Fwd.generation!(rng, pop, idx, ts)
end

pts = Fwd.to_tskit(Fwd.reverse_relabel(ts))

sts = Fwd.simplify(ts, pop.nodes)
print(Fwd.draw_text(ts))
print(Fwd.draw_text(sts))
print(pts.simplify(0:2N-1).draw_text())

xa = Fwd.collect_edges(ts, Fwd.youngest(ts))
xb = Fwd.collect_edges(sts, Fwd.youngest(sts))

pts = Fwd.from_tskit(Fwd.to_tskit(Fwd.reverse_relabel(ts)).simplify())


res = map(1:20) do rep
    @info rep
    ts  = Fwd.init_ts(2N, C) 
    x   = [zeros(Bool, L) for _=1:2N]
    pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    for i=1:20N
        idx = sample(rng, 1:N, 2N)
        pop = Fwd.generation!(rng, pop, idx, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
            @test Fwd.check_edges(ts, pop.nodes)
        end
    end
    rts = Fwd.reverse_relabel(ts)
    pyts = Fwd.to_tskit(rts)
    xs = collect(pyts.breakpoints())
    ts = pyts.diversity(mode="branch", windows=xs) ./ 2
    xs, ts
end

plot(size=(300,200))
map(res) do (xs, ys)
    plot!(xs[1:end-1], first.(ys), color=:black, alpha=0.2)
end
plot!(legend=false)
hline!([2N], lw=2)

_xs = map(x->x[2:end], first.(res))
_ys = last.(res)
xs, ys = Fwd.summarize_wins(_xs, _ys)
plot(xs, vec(mean(ys, dims=1)), linetype=:steppre, color=:black)
hline!([2N], lw=2, legend=false)


# ----------------
N = 100
L = 20
C = 0.1
s = 0.06
h = 0.5
u = s*h/5
xs= collect(range(0.01, 0.02, L))
A = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:2N]

res = map(1:500) do rep
    @info rep
    ts  = Fwd.init_ts(2N, C) 
    x   = [zeros(Bool, L) for _=1:2N]
    pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    for i=1:20N
        ws  = Fwd.eval_fitness(pop)
        idx = sample(rng, 1:N, Weights(ws), 2N)
        pop = Fwd.generation!(rng, pop, idx, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    rts = Fwd.reverse_relabel(ts)
    pyts = Fwd.to_tskit(rts)
    xs = collect(pyts.breakpoints())
    ps = pyts.diversity(mode="branch", windows=xs) ./ 2
    xs, ps
end

_xs = map(x->x[2:end], first.(res))
_ys = last.(res)
xw, yw = Fwd.summarize_wins(_xs, _ys)
plot(xw, vec(mean(yw, dims=1)), linetype=:steppre, color=:black, size=(700,250))
hline!([2N], lw=2, legend=false)
vline!(xs)

# ----------------------
N = 100
L = 1
C = 0.01
s = 0.1
h = 0.5
u = s*h/4
xs= [C/2]
A = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u) for _=1:L], xs)
R = LinearMap(C)
ts= Fwd.init_ts(2N, C) 
x = [zeros(Bool, L) for _=1:2N]
res = map(1:500) do rep
    @info rep
    ts  = Fwd.init_ts(2N, C) 
    x   = [zeros(Bool, L) for _=1:2N]
    pop = Fwd.DiploidWFPopulation(N=N, arch=A, recmap=R, x=x, nodes=collect(1:2N))
    for i=1:20N
        ws  = Fwd.eval_fitness(pop)
        idx = sample(rng, 1:N, Weights(ws), 2N)
        pop = Fwd.generation!(rng, pop, idx, ts)
        if i % 100 == 0
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    rts = Fwd.reverse_relabel(ts)
    pyts = Fwd.to_tskit(rts)
    xs = collect(pyts.breakpoints())
    ps = pyts.diversity(mode="branch", windows=xs) ./ 2
    xs, ps
end

_xs = map(x->x[2:end], first.(res))
_ys = last.(res)
xw, yw = Fwd.summarize_wins(_xs, _ys)
plot(xw, vec(mean(yw, dims=1)), linetype=:steppre, color=:black, size=(700,250))
hline!([2N], lw=2, legend=false)
vline!(xs)

