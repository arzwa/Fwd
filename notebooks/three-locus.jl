
using Distributed
@everywhere using Fwd, TreeSequences, StatsBase
@everywhere using WrightDistribution, Barriers
using Plots; plotsdefault()
using Serialization

const states=[[0,0,0],[0,0,1],[0,1,0],[1,0,0],[0,1,1],[1,0,1],[1,1,0],[1,1,1]]

@everywhere function threelocus_cb(pop)
    pm = proportionmap(pop.x)
    [haskey(pm, x) ? pm[x] : 0.0 for x in states]
end
    
NB  = 500
s   = 10/NB
u   = s/1000
m   = s*0.7
r   = s*0.5
d   = Fwd.distance(r)
L   = 3
flank = L*d
xs  = flank .+ [k*d for k=0:L-1]
C   = last(xs) + flank
R   = LinearMap(C)
AA  = Architecture([BiAllelic(0.0) for _=1:L], xs, R)
AB  = Architecture([BiAllelic(u) for _=1:L], xs, R)
MA  = GPMap([HaploidLocus(0.0, i) for i=1:L])
MB  = GPMap([HaploidLocus(-s, i) for i=1:L])
NA  = 1
nA = collect(1:NA)
nB = collect(1:NB) .+ NA
xA = [ ones(Bool, L) for _=1:NA]
xB = [zeros(Bool, L) for _=1:NB]
popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, x=xA, nodes=nA)
popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, x=xB, nodes=nB)
ngen = 5*10^6
mpop = Fwd.TwoPopOneWay(m, deepcopy(popA), deepcopy(popB))
mpop, ts, qs = simulate!(
    mpop, init_ts(mpop), ngen, x->threelocus_cb(x.popB), every=ngen÷5000)

plot(TS.diffdiv(TS._add_grand_ancestor(ts))[[1,4]], yscale=:log10)

function threelocus_stats(x)
    z = zeros(3)
    for k=1:3
        z[k] = sum(x[filter(i->states[i][k] == 0, 1:length(states))])
    end
    # add LD...
    z 
end

P = pcat(map(threelocus_stats, qs)...)
plot(P)

function g3(p, s, r_1, r_2) 
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_1*p_2^2*s_1*s_2^2 + p_1*p_2*p_3*s_1*s_2*s_3 + p_1*p_2*r_1*s_1*s_2 + p_1*p_2*r_2*s_1*s_2 + p_1*r_1*r_2*s_1 + p_2^3*s_2^3 + p_2^2*p_3*s_2^2*s_3 + 2*p_2^2*r_1*s_2^2 + 2*p_2^2*r_2*s_2^2 + p_2*p_3*r_1*s_2*s_3 + p_2*p_3*r_2*s_2*s_3 + p_2*r_1^2*s_2 + 3*p_2*r_1*r_2*s_2 + p_2*r_2^2*s_2 + p_3*r_1*r_2*s_3 + r_1^2*r_2 + r_1*r_2^2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
end

function g3b(p, s, r_1, r_2)
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_2*s_2 + r_1)*(p_2*s_2 + r_2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_2*s_2 + p_3*s_3 + r_2))
end

function g4(p, s, r_1, r_2)
    p_1, p_2, p_3 = p
    s_1, s_2, s_3 = s
    (p_1*s_1 + r_1)*(p_1*s_1 + p_2*s_2 + r_1 + r_2)/((p_1*s_1 + p_2*s_2 + r_1)*(p_1*s_1 + p_2*s_2 + p_3*s_3 + r_1 + r_2))
end

gs1 = map(p-> g3(p, [s,s,s], r, r), eachrow(P))
gs2 = map(p-> g3b(p, [s,s,s], r, r), eachrow(P))
gs3 = map(p-> g4(p, [s,s,s], r, r), eachrow(P))

scatter(gs1, gs2)
plot!(x->x)

d1 = Wright(2NB*s, NB*(u + mean(gs1)*m), NB*u, 0.5)
d2 = Wright(2NB*s, NB*(u + mean(gs2)*m), NB*u, 0.5)
d3 = Wright(2NB*s, NB*(u + mean(gs3)*m), NB*u, 0.5)
stephist(1 .- P, bins=50, norm=true, label="", legend=:topright)
plot!(range(0,1,200), x->pdf(d1, x), label="derived")
plot!(range(0,1,200), x->pdf(d2, x), label="guess")
plot!(range(0,1,200), x->pdf(d3, x), label="left")


