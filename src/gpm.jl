# Refactor from an `Architecture` point of view to a `GPM` point of view.
# This assumes a quantitative genetic decomposition of the genotypic value.
abstract type Component end

struct GPMap{C<:Component}
    components :: Vector{C}
end
Base.length(m::GPMap) = length(m.components)
Base.getindex(m::GPMap, i) = m.components[i]

GPMap() = GPMap{GenericComponent}(GenericComponent[])

# different components combine additively to yield a phenotype (e.g.
# log-fitness)
function phenotype(gpm::GPMap, x)
    mapreduce(c->eval_component(c, x), +, gpm.components)
end

# a generic component is just a function of genotype
struct GenericComponent <: Component
    fun :: Function
end
eval_component(c::GenericComponent, x) = c.fun(x)

# bi-allelic or integer-valued locus
struct HaploidLocus{T} <: Component
    s :: T
    i :: Int  # index of locus
end
eval_component(c::HaploidLocus, x) = c.s*x[c.i]

# Two locus system
struct HaploidTwoLocus{T} <: Component
    s01 :: T
    s10 :: T
    s11 :: T
    i :: Int
    j :: Int
end
function eval_component(c::HaploidTwoLocus, x)  
    xi, xj = x[c.i], x[c.j]
    return if xi == xj == 0
        0.
    elseif xi == xj == 1
        c.s11
    elseif xi == 1
        c.s10
    else
        c.s01
    end
end 

struct DiploidDominantDMI{T} <: Component
    s :: T
    i :: Int
    j :: Int
end
function eval_component(c::DiploidDominantDMI, x)  
    x1, x2 = x
    cond = (x1[1] == 1 || x2[1] == 1) && (x1[2] == 1 || x2[2] == 1)
    c.s*cond
end

