abstract type MainlandPopulation end

struct HWLEMainland{T,V} <: MainlandPopulation
    m :: T
    p :: V 
end

function randmigrant(rng::AbstractRNG, 
        ml::HWLEMainland{T,V}) where {T,V<:SparseVector}
    SparseVector(ml.p.n, ml.p.nzind, 
        [rand(rng) < ml.p[i] for i in ml.p.nzind])
end

function migration!(rng::AbstractRNG, 
        pop::HaploidWFPopulation, ml::HWLEMainland)
    @unpack m, p = ml
    M = min(pop.N, rand(rng, Poisson(m*pop.N)))
    for i=1:M
        pop.x[i] = randmigrant(rng, ml)
    end
end

function generation!(rng::AbstractRNG, 
        pop::HaploidWFPopulation,
        mainland::MainlandPopulation)
    migration!(rng, pop, mainland)
    generation!(rng, pop)
end
