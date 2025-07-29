struct Edge{T,V}
    parent :: T
    child  :: T
    left   :: V
    rght   :: V
end

function Base.show(io::IO, e::Edge)
    @unpack parent, child, left, rght = e
    write(io, @sprintf("Edge(%d, %d, %.2f, %.2f)", 
        parent, child, left, rght))
end

struct TreeSequence{T<:Integer,V<:Real,W<:Real}
    nodes    :: Vector{W}          # stores birth times 
    edges    :: Vector{Edge{T,V}}  # inheritance edges 
    children :: Vector{Vector{T}}  # for each node, IDs of edges to children
    L        :: V                  # chromosome length
    forward  :: Bool               # forward or backward ts?
end

function youngest(ts::TreeSequence) 
    a, b = extrema(ts.nodes)
    xs = ts.forward ? 
        (findfirst(x->x==b, ts.nodes):length(ts.nodes)) : 
        (1:findlast(x->x==a, ts.nodes))
    collect(xs)
end

function init_ts(N::Int, L::V; forward=true) where V
    nodes = zeros(Int, N)
    edges = Edge{Int,V}[]
    children = [Int[] for _=1:N]
    TreeSequence(nodes, edges, children, L, forward) 
end

function Base.show(io::IO, ts::TreeSequence)
    @unpack nodes, edges, L, forward = ts
    nv = length(nodes)
    ne = length(edges)
    write(io, "TreeSequence(nv=$nv, ne=$ne, L=$L, forward=$forward)")
end

function addnode!(nodes, children::Vector{Vector{T}}, time) where T
    push!(nodes, time)
    push!(children, T[]) 
    return length(nodes)
end

function addnodes!(ts::TreeSequence, N, t)
    x = length(ts.nodes)
    for _=1:N
        addnode!(ts.nodes, ts.children, t)
    end
    return x+1:x+N
end

function addedge!(edges, children, edge)
    push!(edges, edge)
    n = length(edges)
    push!(children[edge.parent], n) 
    return n
end

function addedges!(ts::TreeSequence, edges)
    for e in edges
        addedge!(ts.edges, ts.children, e)
    end
end

function recombine(rng, p1, p2, c, recmap::RecombinationMap)
    bps = rand_breakpoints(rng, recmap)
    edges = map(enumerate(bps)) do (i,x1)
        x0 = (i == 1 ? zero(x1) : bps[i-1])
        isodd(i) ? Edge(p1, c, x0, x1) : Edge(p2, c, x0, x1)
    end
    bps, edges
end

# For simplification algorithm
struct Segment{T,V}
    node :: T
    left :: V
    rght :: V
end

Base.isless(x::Segment, y::Segment) = x.left < y.left

function simplify(ts::TreeSequence{T,V,W}, smpl::Vector{T}) where {T,V,W}
    @unpack nodes, edges, children, L, forward = ts
    #@assert !forward
    nnodes = W[]
    nchildren = Vector{T}[]
    nedges = Edge{T,V}[]
    A = [Segment{T,V}[] for _=1:length(nodes)]
    Q = MutableBinaryMinHeap{Segment{T,V}}()
    for u in smpl
        v = addnode!(nnodes, nchildren, nodes[u])
        A[u] = [Segment(v, 0.0, L)]
    end
    # now nnodes consists of the sample nodes, ordered by their new id's
    # A contains at the index of the old id's the associated active
    # segments (i.e. full genome)
    nn = length(nodes)
    order = forward ? reverse(1:nn) : (1:nn)
    for u in order  # for each node in the original ts, from present to past
        # get the edges of which u is a parent node
        uedges = children[u]  #filter(x->x.parent == u.id, edges)   # S2
        v = -1
        for k in uedges  # S3
            e = edges[k]
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
                    v = addnode!(nnodes, nchildren, nodes[u])
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
            push!(A[u], α)
        end
    end
    sort!(nedges, by=x->(x.parent, x.child, x.rght, x.left))
    nnedges = Edge{T,V}[]
    start = 1
    for j in 2:length(nedges)
        ei = nedges[j-1]
        ej = nedges[j]
        condition = (
            ei.rght != ej.left || ei.parent != ej.parent || ei.child != ej.child)
        if condition  # edges must not be merged
            ek = Edge(ei.parent, ei.child, nedges[start].left, ei.rght)
            addedge!(nnedges, nchildren, ek)
            start = j
        end
    end
    if length(nedges) > 0 
        ei = last(nedges)
        ek = Edge(ei.parent, ei.child, nedges[start].left, ei.rght)
        addedge!(nnedges, nchildren, ek)
    end
    if forward
        maxt = maximum(nnodes)
        sts = reverse_relabel(TreeSequence(nnodes, nnedges, nchildren, L, false))
        sts.nodes .= maxt .- sts.nodes
    else
        sts = TreeSequence(nnodes, nnedges, nchildren, L, false)
    end
    return sts
end

function to_tskit(ts::TreeSequence)
    @unpack nodes, edges, L, forward = ts
    @assert !forward
    tbc = tskit.TableCollection(L)
    for t in nodes
        tbc.nodes.add_row(time=t, flags= t==0 ? 1 : 0)
    end
    for e in edges
        tbc.edges.add_row(e.left, e.rght, e.parent-1, e.child-1)
        # don't forget: 0-based indexing in Python!
    end
    tbc.sort()
    tbc.tree_sequence()
end

function from_tskit(ts)
    tabs = ts.tables
    nodes = map(x->x.time, tabs.nodes)
    edges = map(x->Edge(x.parent+1, x.child+1, x.left, x.right), tabs.edges)
    children = [filter(i->edges[i].parent == j, 1:length(edges)) for j=1:length(nodes)]
    L = ts.sequence_length
    TreeSequence(nodes, edges, children, L, false)
end

# nodes stay the same, edges do not change internally either, only the edge
# index and hence `children`
function sort_edges(ts::TreeSequence)
    e, c = sort_edges(ts.edges, ts.children)
    reconstruct(ts, edges=e, children=c)
end

function sort_edges(edges, children)
    o = sortperm(edges, by=x->(x.parent, x.child, x.left, x.rght))
    i = invperm(o)
    nchildren = map(children) do xs
        [i[x] for x in xs]
    end
    edges[o], nchildren 
end

function reverse_relabel(ts::TreeSequence)
    # assumes currently sorted from past to present, and that present is
    # maximum T, will yield the reverse: i.e. sorted from present to past,
    # with the present being t=0
    @unpack nodes, edges, children, L, forward = ts
    T = forward ? nodes[end] : nodes[1]
    n = length(nodes)
    nedges = map(edges) do e
        Edge(n - (e.parent-1), n - (e.child-1), e.left, e.rght)
    end
    e, c = sort_edges(nedges, children)
    TreeSequence(reverse(T .- nodes), e, reverse(c), L, !forward)
end

function collect_edges(ts, nodes)
    es = [filter(x->x.child == n, ts.edges) for n in nodes]
end

# XXX note that not all leaf nodes will have edges for their full sequence
# length in a simplified ts when the sample has not fully coalesced...
function check_edges(ts, nodes)
    es = collect_edges(ts, nodes)
    ls = map(ex->sum([x.rght-x.left for x in ex]), es)
    all(ls .== ts.L)
end

function check_children(ts)
    @unpack children, edges = ts
    map(1:length(children)) do i
        all([edges[k].parent == i for k in children[i]])
    end |> all
end

draw_text(ts) = ts.draw_text()
function draw_text(ts::TreeSequence)
    if ts.forward 
        ts = reverse_relabel(ts)
    end
    to_tskit(ts).draw_text()
end

# neutral WF -----------------------------------------------------------
# pop is a collection of node IDs, where index k and N+k give the two
# haplotypes in individual k
function generation!(rng, pop, ts::TreeSequence{T}, recmap) where T
    @unpack nodes, edges, children, L = ts
    n = length(ts.nodes) 
    N = length(pop)÷2
    t = nodes[pop[1]] + 1
    pop′ = [n+i for i=1:2N]
    for k=1:N
        # choose parent for first haplotype
        i = rand(rng, 1:N)
        eik = _recombine(rng, pop[i], pop[N+i], pop′[k], recmap) 
        for e in eik
            addedge!(edges, children, e)
        end
        # choose parent for the second haplotype
        j = rand(rng, 1:N)
        ejk = _recombine(rng, pop[j], pop[N+j], pop′[N+k], recmap) 
        for e in ejk
            addedge!(edges, children, e)
        end
    end
    nodes = vcat(nodes, fill(t,2N))
    children = vcat(children, [T[] for _=1:2N])
    pop′, TreeSequence(nodes, edges, children, L, true)
end

function _recombine(rng, p1, p2, c, recmap::RecombinationMap)
    bps = rand_breakpoints(rng, recmap)
    map(enumerate(bps)) do (i,x1)
        x0 = i == 1 ? zero(x1) : bps[i-1]
        isodd(i) ? Edge(p1, c, x0, x1) : Edge(p2, c, x0, x1)
    end
end

