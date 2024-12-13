using Pkg; Pkg.activate("/home/arthur_z/dev/Fwd")
using Test, Fwd, BenchmarkTools, Profile, ProfileView
using WrightDistribution, Barriers, Windows
using Plots; plotsdefault()
using Fwd: kb, Mb
rng = default_rng()

@testset "" begin
    L = 1kb
    x1 = spzeros(Bool, L)
    x2 = spzeros(Bool, L)
    x2[501:end] .= true
    x3 = Fwd.recombination(x1, x2, [250, 500, 750])
    @test x3[1:250]   == x1[1:250]
    @test x3[251:500] == x2[251:500]
    @test x3[501:750] == x1[501:750]
    @test x3[751:end] == x2[751:end]
end

# A single locus
loci = [Fwd.HaploidLocus(-0.005, 0.002)]
xs   = [Fwd.Region([1])]
arch = Fwd.Architecture(loci, xs)
N    = 200
X    = [spzeros(Bool, 1) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(0.1))
ps   = Fwd.simulate(rng, pop, 5000N)

Y = vcat(map(first, ps)...)
d = Wright(-N*(2loci[1].s), N*loci[1].u, N*loci[1].u, 0.5)
P = stephist(Y, bins=0:0.025:1, norm=true)
pxs = 0:0.025:1
plot!(pxs, map(p->pdf(d, p), pxs))

# Two loci
loci = [Fwd.HaploidLocus(-0.005, 0.002), Fwd.HaploidLocus(0.0, 0.002)]
xs   = [Fwd.Region([1]), Fwd.Region([2])]
arch = Fwd.Architecture(loci, xs)
N    = 200
X    = [spzeros(Bool, 2) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(0.1))
ps   = Fwd.simulate(rng, pop, 1000N)

Y1 = vcat(map(first, ps)...)
Y2 = vcat(map( last, ps)...)
d = Wright(-N*(2loci[1].s), N*loci[1].u, N*loci[1].u, 0.5)
P1 = stephist(Y1, bins=0:0.025:1, norm=true)
pxs = 0:0.025:1
plot!(pxs, map(p->pdf(d, p), pxs))
d = Wright(-N*(2loci[2].s), N*loci[2].u, N*loci[2].u, 0.5)
P2 = stephist(Y2, bins=0:0.025:1, norm=true)
pxs = 0:0.025:1
plot!(pxs, map(p->pdf(d, p), pxs))
plot(P1, P2, size=(600,200))


# bunch of loci, equal effect
L    = 1kb
N    = 500
loci = [Fwd.HaploidLocus(-0.005, 2/L)]
xs   = [Fwd.Region(1:10:1kb)]
arch = Fwd.Architecture(loci, xs)
X    = [spzeros(Bool, L) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(0.0))
ps   = Fwd.simulate(rng, pop, 100N)

@time Fwd.generation!(rng, pop);

Y = hcat(map(p->[p[i] for i in arch.xs[1].xs], ps)...)
d = Wright(-2N*loci[1].s, 2N*loci[1].u, 2N*loci[1].u, 0.5)
P = stephist(vcat(Y[:,10N:end]...), bins=0:0.033:1, norm=true,
    color=:lightgray, fill=true, label="simulation")
pxs = 0:0.001:1
plot!(pxs, map(p->pdf(d, p), pxs), color=:black, ylabel="density",
    xlabel="\$q\$", label="single locus prediction")
plot!(title="\$r=$(pop.recmap.r), s=$(loci[1].s), u=$(loci[1].u)\$", size=(340,250))

# selected and neutral in between (BGS model)
L    = 1kb
ls   = 451:550
ln   = [1:450, 551:1kb]
N    = 500
Ud   = 1.0
Un   = 0.001
loci = [
    Fwd.HaploidLocus(-0.02, Ud/length(ls)), 
    Fwd.HaploidLocus(0.0, Un/length(ln))]
xs   = [Fwd.Region(ls), Fwd.Region(ln)]
arch = Fwd.Architecture(loci, xs)
X    = [spzeros(Bool, L) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(0.1/L))
ps   = Fwd.simulate(rng, pop, 100N)

Yn = hcat(map(p->Array([p[i] for i in arch.xs[2].xs][1]), ps)...)
Ys = hcat(map(p->Array([p[i] for i in arch.xs[1].xs]), ps)...)

bins = 0:0.04:1
stephist(vcat(Yn[1:10,10N:end]...), bins=bins)
stephist!(vcat(Yn[end-10:end,10N:end]...), bins=bins)
stephist!(vcat(Ys[:,10N:end]...), bins=bins)

stephist(vcat(Yn[1:10,10N:end]...), bins=0:0.025:1, norm=true)
stephist!(vcat(Yn[end-10:end,10N:end]...), bins=0:0.025:1, norm=true)
pxs = 0:0.025:1
d = Wright(0., 2N*loci[2].u, 2N*loci[2].u, 0.5)
plot!(pxs, map(p->pdf(d, p), pxs))


# Large neutral model
L    = 5Mb
N    = 500
u    = 1e-9
r    = 1e-8
loci = [Fwd.HaploidLocus(0.0, u)]
xs   = [Fwd.Region(1:L)]
arch = Fwd.Architecture(loci, xs)
X    = [spzeros(Bool, L) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(r))
ps   = Fwd.simulate(rng, pop, 100N, remove_fixed=100)

@btime Fwd.simulate(rng, pop, 10N, remove_fixed=100, cb=(x)->NaN, show_progress=false);


# Mainland island migration
Ls   = 1.0
L    = 100
s    = Ls/L
Ns   = 10.0
N    = ceil(Int, Ns/s)
u    = s/200
m    = s*0.5
r    = 0.1
loci = [Fwd.HaploidLocus(-s, u)]
xs   = [Fwd.Region(1:L)]
arch = Fwd.Architecture(loci, xs)
X    = [spzeros(Bool, L) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(r))
ml   = Fwd.HWLEMainland(m, sparse(ones(L)))
ps   = Fwd.simulate(rng, pop, 10N, ml)
p̄    = 1 .- mean(ps)
mean(p̄)
# should be ~ 0.93 (checked with Sewall)

using Barriers
recmap= Barriers.linearmap(L, Fwd.invhaldane(r)*(L-1))
ys    = [recmap[i] for i=1:L]
locus = Barriers.MILocus(2s, 0.5, u, 0.0, float(N))
Ma    = Barriers.MIModel(ml.m, ys, [locus for _=1:L])
mean(Ma.Ep)
plot(p̄)
plot!(Ma.Ep)


# --------------------------------------------------------
# BGS + divergent selection
G = 1Mb
r = 1e-6

## BGS
U  = 0.2
sd = 0.005
Gd = 500kb
Ld = 1000
ud = U/Ld 
xd = Fwd.Region(((G - Gd)÷2 + 1):(Gd÷Ld):((G + Gd)÷2))
ld = Fwd.HaploidLocus(-sd, ud)

## Divergent seln.
L  = 20 
Ls = 0.2
sa = Ls/L
ua = sa/200
Ns = 20.0
N  = ceil(Int, Ns/sa)
Δ  = G÷L
xa = Fwd.Region(Δ÷2:Δ:G)
la = Fwd.HaploidLocus(-sa, ua)
m  = sa*0.7
p̄  = spzeros(G)
p̄[xa.xs] .= 1.0
p̄[xd.xs] .= ud/sd
ml = Fwd.HWLEMainland(m, p̄)
ys = [recmap[i] for i in xa.xs]
@info (ys[2]-ys[1])/sa

## Neutral alleles
Nn = 1000
xn = setdiff(1:G, union(xd.xs, xa.xs))
xn = Fwd.Region(ceil.(Int, range(1, length(xn), Nn)))
un = 1/xn.L
ln = Fwd.HaploidLocus(0.0, un)
p̄[xn.xs] .= 0.5
ml = Fwd.HWLEMainland(m, p̄)

## theoretical prediction
using Barriers
recmap= Barriers.linearmap(G, Fwd.invhaldane(r)*(G-1))
locus = Barriers.MILocus(2sa, 0.5, ua, 0.0, float(N))
Ma    = Barriers.MIModel(ml.m, ys, [locus for _=1:L])
yd    = [recmap[i] for i in xd.xs]
Md    = Barriers.BGSModel(yd, [Barriers.DelLocus(sd, ud) for _=1:Ld])
Bs    = map(x->Barriers.B( Md, recmap[x]), xn.xs)
mes   = map(x->Barriers.me(Ma, recmap[x]), xn.xs)
fstpred = 1 ./ (1 .+ 2N .* Bs .* mes)
plot(fstpred, size=(500,200), color=:salmon, ylabel="\$F_{\\mathrm{ST}}\$",
    fill=true, fillalpha=0.2)

arch = Fwd.Architecture([ld, la, ln], [xd, xa, xn])
X    = [spzeros(Bool, G) for i=1:N]
pop  = Fwd.HaploidWFPopulation(X, arch, Fwd.LinearMap(r))

ps   = Fwd.simulate(rng, pop, 10N, ml)

plot(fstpred, size=(500,200), color=:salmon, ylabel="\$F_{\\mathrm{ST}}\$",
    fill=true, fillalpha=0.2)
pns = map(p->p[xn.xs], ps[2N:10:end])
fst = mean(map(p-> 1 .- (p .* (1 .- p)) ./ 0.25, pns))
plot!(midwindow(mean, fst, 10, 5)..., color=:black, size=(500,200))

plot(midwindow(mean, fst, 2, 1)..., color=:black, size=(500,200))
plot!(fstpred, size=(500,200), color=:salmon, ylabel="\$F_{\\mathrm{ST}}\$")

