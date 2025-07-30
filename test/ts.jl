using Fwd, Test, PyCall
import Fwd: Node, Edge, TreeSequence, LinearMap
using Random
rng = Random.default_rng()


Random.seed!(12)
N = 1000
pop = [k for k=1:2N]
nodes = [Fwd.Node(0.0, 1) for i=1:2N]
edges = Edge{Int64,Float64}[]
children = [Int64[] for _=1:2N]
L = 0.1
recmap = LinearMap(L)
ts = TreeSequence(nodes, edges, children, L, true)
ngen = 10
for _=1:ngen
    pop, ts = _generation!(rng, pop, ts, recmap)
end
sts = Fwd.simplify(ts, collect((ngen+1)*2N-10:(ngen+1)*2N))

function _generation!(rng, pop, ts::TreeSequence{T}, recmap) where T
    @unpack nodes, edges, children, L = ts
    n = length(ts.nodes) 
    N = length(pop)÷2
    t = nodes[pop[1]].time + 1
    pop′ = [n+i for i=1:2N]
    for k=1:N
        # choose parent for first haplotype
        i = rand(rng, 1:N)
        eik = Fwd._recombine(rng, pop[i], pop[N+i], pop′[k], recmap) 
        for e in eik
            Fwd.addedge!(edges, children, e)
        end
        # choose parent for the second haplotype
        j = rand(rng, 1:N)
        ejk = Fwd._recombine(rng, pop[j], pop[N+j], pop′[N+k], recmap) 
        for e in ejk
            Fwd.addedge!(edges, children, e)
        end
    end
    nodes = vcat(nodes, [Fwd.Node(t, 1) for _=1:2N])
    children = vcat(children, [T[] for _=1:2N])
    pop′, TreeSequence(nodes, edges, children, L, true)
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

