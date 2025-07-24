# Two (explicit) populations model, with unidirectional migration and **no
# selection**. This is to check our neutral predictions for the case with
# finite populations and unidirectional migration (IM-without-I model).
using Fwd, Plots, Test, BenchmarkTools, WrightDistribution, Distributions
using Serialization, Barriers
default(grid=false, fontfamily="Computer modern", framestyle=:box,
    legend=false, fg_legend=:transparent)

# Parameters --------------------------------------------------------
NA  = 500         # population size
NB  = 100
m   = 5e-3
u   = m/100
L   = 1000

fstfun(m, NA, NB)

# Define the population model ---------------------------------------
# Seed the RNG
rng     = Random.seed!(99)
arch    = Architecture([HaploidBiLocus(0.0, u) for i=1:L], fill(Inf, L))
hapsA   = map(_->Vector{Bool}(rand(L) .< 0.5), 1:NA)
hapsB   = map(_->Vector{Bool}(rand(L) .< 0.5), 1:NB)
recmap  = Fwd.Unlinked()
island  = HaploidWFPopulation(hapsB, arch, recmap)
mainland= HaploidWFPopulation(hapsA, arch, recmap) 
metapop = MainlandIsland(island, mainland, m)

# Conduct a simulation ----------------------------------------------
Xs = simulate!(metapop, 200NA, every=20, 
    callback=P->(mean(P.island.x), mean(P.mainland.x)))

# Some plots
plot(mean(first.(Xs) .- last.(Xs)), ylabel="allelic divergence")

Xi = hcat(first.(Xs)...)
Xm = hcat(last.(Xs)...)

plot()
map(rand(1:L, 1)) do i
    plot!(Xi[i,:]); plot!(Xm[i,:])
end; plot!(size=(800,200))

map(rand(1:L, 9)) do i
    stephist( Xm[i,:], bins=0:0.01:1.01)
    stephist!(Xi[i,:], bins=0:0.01:1.01)
end |> x->plot(x...)

qi = Xi[:,end]
qm = Xm[:,end]


fim = map(Xs) do (qi, qm)
    pii = mean(qi .* (1 .- qi))
    pim = mean(qm .* (1 .- qm))
    pib = mean(qi .* (1 .- qm) + (1 .- qi) .* qm)/2
    piw = (pii + pim)/2
    (pib - piw)/(pib + piw)
end 

stephist(fim)
vline!([fstfun(m, NA, NB)])
vline!([mean(fim)])

plot(fim)
hline!([fstfun(m, NA, NB)])
hline!([mean(fim)])

fstfun = eval(
    :(function (m, N_A, N_B)
      #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =# @inbounds begin
              #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =#
              begin
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:389 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:390 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:391 =#
                  (/)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(1//2, N_A), m)), (*)((*)(1//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-1//1, N_A), (^)(m, 2))), (*)((*)(-1//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)(1//2, N_A), (^)(m, 3))), (*)((*)(1//2, N_B), (^)(m, 3))), (+)((+)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(3//2, N_A), m)), (*)((*)(7//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-3//1, N_A), (^)(m, 2))), (*)((*)(-5//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)((*)(4//1, N_A), N_B), (^)(m, 2))), (*)((*)(3//2, N_A), (^)(m, 3))), (*)((*)(3//2, N_B), (^)(m, 3))), (*)((*)((*)(-2//1, N_A), N_B), (^)(m, 3))))
              end
          end
  end)
)

