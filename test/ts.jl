using Fwd, Test, PyCall
import Fwd: Node, Edge, TreeSequence, LinearMap
using Random
rng = Random.default_rng()

Random.seed!(12)
N = 1000
pop = [k for k=1:2N]
nodes = fill(0, 2N)
edges = Edge{Int64,Float64}[]
children = [Int64[] for _=1:2N]
L = 0.1
recmap = LinearMap(L)
ts = TreeSequence(nodes, edges, children, L, true)
ngen = 1000
for _=1:ngen
    pop, ts = Fwd.generation!(rng, pop, ts, recmap)
end

rts = Fwd.reverse_relabel(ts)
sts = Fwd.simplify(rts, findall(iszero, rts.nodes))

fts = Fwd.simplify(ts, pop)

sts = Fwd.sort_edges(ts)
rts = Fwd.reverse_relabel(sts)
rrts = Fwd.reverse_relabel(rts)
@test all(sts.edges .== rrts.edges)
@test all(sts.children .== rrts.children)

smpl = [i for i=1:length(ts.nodes) if ts.nodes[i] == ngen]
ssts = Fwd.simplify(sts, smpl)

smpl = [i for i=1:length(sts.nodes) if sts.nodes[i] == 0][1:5]
ssts = Fwd.simplify(rts, smpl)
tspy = Fwd.to_tskit(ssts)
tspy.num_trees
print(tspy.at(0.05).draw_text())

