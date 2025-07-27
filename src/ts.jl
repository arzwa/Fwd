
struct Node{T}
    id   :: T
    time :: Int64
end

struct Edge{T,V}
    parent :: T
    child  :: T
    left   :: V
    rght   :: V
end

function Base.show(io::IO, e::Edge)
    s = @sprintf "Edge(%d, %d, %.2f, %.2f)" e.parent e.child e.left e.rght
    write(io, s)
end

struct TreeSequence{T,V}
    nodes :: Vector{Node{T}}
    edges :: Vector{Edge{T,V}}
    L :: V 
end

function addnode!(nodes, time)
    n = length(nodes)
    push!(nodes, Node(n+1, time))
    n + 1
end

Diploid{T} = Tuple{Node{T},Node{T}}
Haploid{T} = Node{T}

struct Segment{T,V}
    node :: T
    left :: V
    rght :: V
end

Base.isless(x::Segment, y::Segment) = x.left < y.left

function simplify(ts::TreeSequence{T,V}, smpl::Vector{Node{T}}) where {T,V}
    @unpack nodes, edges, L = ts
    nnodes = Node{T}[]
    nedges = Edge{T,V}[]
    A = [Segment{T,V}[] for _=1:length(nodes)]
    Q = MutableBinaryMinHeap{Segment{T,V}}()
    for u in smpl
        v = addnode!(nnodes, u.time)
        A[u.id] = [Segment(v, 0.0, L)]
    end
    # now nnodes consists of the sample nodes, ordered by their new id's
    # A contains at the index of the old id's the associated active
    # segments (i.e. full genome)
    for u in sort(nodes, by=x->x.time)  # for each node in the original ts XXX sorted from present to past
        # get the edges of which u is a parent node
        # XXX this would be much faster if one had an actual graph data
        # structure?
        uedges = filter(x->x.parent == u.id, edges)   # S2
        v = -1
        for e in uedges  # S3
            for x in A[e.child]
                if x.rght > e.left && e.rght > x.left
                    # we found a segment of e.child overlapping with the
                    # edge (i.e. what e.child inherited from u)
                    y = Segment(x.node, max(x.left, e.left), min(x.rght, e.rght))
                    push!(Q, y)
                end
            end
        end
        #@info "before while (S4)" u uedges filter(!isempty, A) Q
        while !isempty(Q)  # S4
            l = first(Q).left
            r = L 
            X = Segment{T,V}[]
            while !isempty(Q) && first(Q).left == l
                x = pop!(Q)
                push!(X, x)
                r = min(r, x.rght)
            end
            r = !isempty(Q) ? min(r, first(Q).left) : r
            if length(X) == 1  # S5, no overlapping segments
                x = X[1]
                α = x
                if !isempty(Q) && first(Q).left < x.rght
                    #α = Segment(x.node, x.left, first(Q).rght)
                    α = Segment(x.node, x.left, first(Q).left)
                    x = Segment(x.node, first(Q).left, x.rght)
                    push!(Q, x)
                end
            else  #  S6, there's overlap, new output node
                if v == -1
                    v = addnode!(nnodes, u.time)
                end
                # S7 record edges
                α = Segment(v, l, r)
                for x in X
                    push!(nedges, Edge(v, x.node, l, r))
                    if x.rght > r
                        x = Segment(x.node, r, x.rght)
                        push!(Q, x)
                    end
                end
            end
            push!(A[u.id], α)
        end
    end
    sort!(nedges, by=x->(x.parent, x.child, x.rght, x.left))
    nnedges = Edge{T,V}[]
    start = 1
    for j in 2:length(nedges)
        ei = nedges[j-1]
        ej = nedges[j]
        if (ei.rght != ej.left || 
                ei.parent != ej.parent || 
                ei.child != ej.child)
            ek = Edge(ei.parent, ei.child, nedges[start].left, ei.rght)
            push!(nnedges, ek)
            start = j
        end
    end
    ei = last(nedges)
    ek = Edge(ei.parent, ei.child, nedges[start].left, ei.rght)
    TreeSequence(nnodes, nnedges, L)
end

function recombine(rng, p1, p2, c, recmap::RecombinationMap)
    bps = rand_breakpoints(rng, recmap)
    map(enumerate(bps)) do (i,x1)
        x0 = i == 1 ? zero(x1) : bps[i-1]
        isodd(i) ? Edge(p1.id, c.id, x0, x1) : Edge(p2.id, c.id, x0, x1)
    end
end

# neutral WF
function generation!(rng, pop, ts::TreeSequence, recmap)
    @unpack nodes, edges, L = ts
    k = length(ts.nodes) 
    N = length(pop)
    t = pop[1][1].time + 1
    pop′ = [(Node(k+i,t), Node(N+k+i,t)) for i=1:N]
    for k=1:N
        # choose parent for first haplotype
        i = rand(rng, 1:N)
        eik = recombine(rng, pop[i]..., pop′[k][1], recmap) 
        # choose parent for the second haplotype
        j = rand(rng, 1:N)
        ejk = recombine(rng, pop[j]..., pop′[k][2], recmap) 
        edges = vcat(edges, eik, ejk)
    end
    nodes = vcat(nodes, first.(pop′), last.(pop′))
    pop′, TreeSequence(nodes, edges, L)
end

function to_tskit(ts::TreeSequence)
    @unpack nodes, edges, L = ts
    tbc = tskit.TableCollection(L)
    sort!(edges, by=x->(x.parent, x.child))
    for n in nodes
        tbc.nodes.add_row(time=n.time, flags=n.time==0 ? 1 : 0)
    end
    for e in edges
        tbc.edges.add_row(e.left, e.rght, e.parent-1, e.child-1)
        # don't forget: 0-based indexing in Python!
    end
    tbc.sort()
    tbc.tree_sequence()
end

function reverse_time(ts::TreeSequence)
    @unpack nodes, edges, L = ts
    T = maximum(x->x.time, nodes)
    n = [Node(n.id, T-n.time) for n in nodes]
    TreeSequence(n, edges, L)
end

function reverse_relabel(ts::TreeSequence)
    # assumes currently sorted from past to present, and that present is
    # maximum T, will yield the reverse: i.e. sorted from present to past,
    # with the present being t=0
    @unpack nodes, edges, L = ts
    T = nodes[end].time
    n = length(nodes)
    nnodes = [Node(k, T-node.time) for (k,node) in enumerate(reverse(nodes))]
    nedges = map(edges) do e
        Edge(n - (e.parent-1), n - (e.child-1), e.left, e.rght)
    end
    sort!(nedges, by=x->(x.parent, x.child, x.left, x.rght))
    TreeSequence(nnodes, nedges, L)
end


