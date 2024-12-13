using Plots
plotsdefault()

annot = "/home/arthur_z/dev/Fwd/data/MpTak1_v7.1.gff"

# extract a chromosome
flter(x) = x[1] == "chr8" && x[3] == "gene"
genes = open(annot, "r") do f
    readlines(f) .|> split |> x->filter(flter, x) .|> x->parse.(Int64, getindex(x, [4,5]))
end
L = 24713864  # chr 8 size

@info issorted(genes)
@info all([genes[i+1][1] > genes[i][2] for i=1:length(genes)-1])
# There are overlapping genes.

# Join overlapping genes
genes_ = Vector{Int64}[genes[1]]
for i=2:length(genes)
    if genes[i][1] < genes[i-1][2]
        genes_[end][2] = genes[i][2]
    else 
        push!(genes_, genes[i])
    end
end
# Get a subset
genes = genes_[1:20]

# Set the length of the cut-off bit of chromosome
L = genes[end][2] + genes[1][1]

plot(size=(1800,200), legend=false)
map(genes) do g
    plot!(g, [1,1], fill=true, color=:black, lw=0)
end
plot!(xlim=(0,L), margin=5Plots.mm)
# This seems managable

# Let's do a theoretical prediction
using Fwd, Barriers, Random
using Fwd: kb, Mb
rng = default_rng()

coding = Fwd.Region([g[1]:3:g[2] for g in genes])
coding.L/L
r     = 1e-6
Ne    = 1000
Ns    = 20.0
ud    = 0.1/coding.L
s     = Ns/Ne
rmap  = Barriers.linearmap(L, r*L)
xd    = union(collect.(coding.xs)...)
yd    = [rmap[i] for i in xd]
Md    = Barriers.BGSModel(yd, [Barriers.DelLocus(s, ud) for _=1:length(yd)])
bx    = 1:1kb:L
Bs    = map(x->Barriers.B( Md, rmap[x]), bx)
plot(size=(900,200), legend=false)
map(genes) do g
    plot!(g, [1,1], fill=true, color=:black, lw=0, alpha=0.3)
end
plot!(xlim=(0,L), margin=5Plots.mm)
plot!(twinx(), bx, Bs, legend=false, ylabel="\$B\$", color=:black)

noncoding = [genes[i][2]+1:genes[i+1][1]-1 for i=1:length(genes)-1]
noncoding = [[1:genes[1][1]-1] ; noncoding; [genes[end][2]+1:L]]
noncoding = [noncoding ; [g[1]+1:3:g[2]-1 for g in genes]]
noncoding = Fwd.Region(noncoding)
noncoding.L
un = 1/noncoding.L
θ  = 2un*Ne

loci = [
    Fwd.HaploidLocus(-s, ud), 
    Fwd.HaploidLocus(0.0, un)
]
arch = Fwd.Architecture(loci, [coding, noncoding])
X    = [spzeros(Bool, L) for i=1:Ne]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(r))
ps   = Fwd.simulate(rng, pop, 5Ne, remove_fixed=100)

ps   = [ps; Fwd.simulate(rng, pop, 5Ne, remove_fixed=100)]

pqs = map(p->p .* (1 .- p), ps[1Ne:end])
pqq = mean(pqs)
#pqq = vcat([pqq[r] for r in noncoding.xs]...)

pqq = Vector{Union{Float64,Missing}}(pqq)
for r in coding.xs
    pqq[r] .= missing
end

plot(size=(600,200), margin=5Plots.mm, legend=:topleft)
#plot!(twinx(), bx, Bs, label="\$B\$", legend=:bottomright, ylabel="\$B\$", color=:black)
map(genes) do g
    plot!(twinx(), g, [1,1], fill=true, color=:black, lw=0, alpha=0.3, legend=false, yticks=false)
end
xx, yy = midwindow(mean, pqq, 200, 200)
plot!(xx, 2yy, label="\$\\pi\$", color=:black, lw=1)
plot!(bx, (2Ne .* Bs * un) ./ (4Ne .* Bs * un .+ 1), label="\$\\mathbb{E}[\\pi]\$", lw=2)
coding = Fwd.Region([g[1]:3:g[2] for g in genes])
hline!([2Ne*un*mean(Bs)], label="\$2N_eu\\overline{B}\$", lw=2)
plot!(size=(900,200))

