# NOTE: everything is written so that any population type could be used
# given that it implements the necessary functions for (1) ts
# initialisation, (2) migration, (3) generation.
"""
    TwoPop

Migration is from A to B forward in time. Migration replaces an expected
proportion `m` of the B population by A individuals, without affecting A.
"""
struct TwoPopOneWay{P,T}
    m    :: T
    popA :: P
    popB :: P
end

active_nodes(m::TwoPopOneWay) = [m.popA.nodes ; m.popB.nodes]

"""
    migration!(rng, metapop)

!!! This migration function implements migration-as-copying, i.e. a
proportion `m` of the B population is replaced by indviduals from A, but A
is unaffected.  (biologically this could correspond to sending out asexual
propagules -- although this is more relevant for haplodiplontic cryptogams
etc. than for diploids).
"""
function migration!(rng, metapop::TwoPopOneWay) 
    @unpack popA, popB, m = metapop
    NB = popB.N
    NA = popA.N
    nmig = min(NB, rand(rng, Poisson(m*NB)))
    idx = sample(rng, 1:NA, nmig) 
    for k=1:nmig
        migrate!(popA, popB, idx[k], k)
    end
end
# Other possibilities for migration:
# 1. as in SLiM, a proportion m of the offspring is replaced by offspring
# from a pair of individuals from the source population
# 2. ...

function generation!(rng, metapop::TwoPopOneWay, ts::TreeSequence)
    migration!(rng, metapop)
    @unpack popA, popB = metapop
    popA_ = generation!(rng, popA, ts, popid=1)
    popB_ = generation!(rng, popB, ts, popid=2)
    reconstruct(metapop, popA=popA_, popB=popB_)
end

function simplify!(pop::TwoPopOneWay, ts::TreeSequence)
    @unpack popA, popB = pop
    ns = active_nodes(pop)
    sts = simplify(ts, ns, keep_roots=true)
    #T = time(sts.nodes[end])
    # XXX the indices are not what I expected?
    #popA.nodes .= findall(x->time(x) == T && population(x) == 1, sts.nodes)
    #popB.nodes .= findall(x->time(x) == T && population(x) == 2, sts.nodes)
    nv = length(sts.nodes)
    NA = length(popA.nodes)
    NB = length(popB.nodes)
    na = (nv-NA-NB+1):(nv-NB)
    nb = (nv-NB+1):nv
    popA.nodes .= collect(na)
    popB.nodes .= collect(nb)
    return pop, sts
end

function init_ts(pop::TwoPopOneWay, L)
    @unpack popA, popB = pop
    tsa = init_ts(popA, L, popid=1)
    tsb = init_ts(popB, L, popid=2)
    ns = [tsa.nodes; tsb.nodes]
    es = tsa.edges
    cs = [tsa.children; tsb.children] 
    TreeSequence(ns, es, cs, L, true)
end

