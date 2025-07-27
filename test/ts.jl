using Fwd
import Fwd: Node, Edge, TreeSequence, LinearMap
using Random
rng = Random.default_rng()

Random.seed!(12)
N = 5
pop = [(Node(k,0), Node(N+k, 0)) for k=1:N]
nodes = vcat(first.(pop), last.(pop))
edges = Edge{Int64,Float64}[]
ts = TreeSequence(nodes, edges, 1.0)
recmap = LinearMap(1.0)
for _=1:8
    pop, ts = Fwd.generation!(rng, pop, ts, recmap)
end
ts = Fwd.reverse_relabel(ts)

sts = Fwd.simplify(ts, [x for x in ts.nodes if x.time == 0][5:8])

tspy = Fwd.to_tskit(sts)

