using Fwd
using Random
rng = Random.default_rng()

N = 1000
s = 0.02
h = 0.5
u = s*h/5
L = 50
C = 0.1
xs= C/2L:C/L:C 
A = Architecture([Fwd.DiploidBiLocus(-s*h, -s, u) for _=1:L], xs)
R = LinearMap(C)
x = [zeros(Bool, L) for _=1:2N]

pop = Fwd.DiploidWFPopulation(N, A, R, x, deepcopy(x))
ts  = Fwd.init_ts(2N, C) 
tsrec = Fwd.TreeSequenceRecorder(collect(1:2N), ts)
for i=1:10000
    pop = Fwd.generation!(rng, pop, tsrec)
    if i % 20 == 0
        ts = Fwd.simplify(tsrec.ts, tsrec.active)
        active = findall(x->x==i, ts.nodes)
        tsrec = Fwd.TreeSequenceRecorder(active, ts)
    end
end
rts = Fwd.reverse_relabel(tsrec.ts)
pyts = Fwd.to_tskit(rts)

xxs = pyts.breakpoints() |> collect
div = pyts.diversity(collect(pyts.samples()), mode="branch", windows=xxs)
mdiv = sum((xxs[2:end] .- xxs[1:end-1]) .* div ./ C)
plot(xxs[2:end], div, linetype=:steppre, color=:black, size=(800,200), legend=false)
hline!([2N, mdiv])
vline!(xs)


