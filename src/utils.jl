# Simulation routines...
simulate!(pop::AbstractPop, args...) = simulate!(default_rng(), args...)
function simulate!(rng::AbstractRNG, pop, ts, ngen; simplify=100)
    @showprogress for t=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        if t % simplify == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    return pop, ts
end

function simulate!(rng::AbstractRNG, pop, ts, ngen, cb::Function; simplify=100)
    ys = [cb(pop)]
    @showprogress for t=1:ngen
        pop = Fwd.generation!(rng, pop, ts);
        push!(ys, cb(pop))
        if t % simplify == 0 
            pop, ts = Fwd.simplify!(pop, ts)
        end
    end
    return pop, ts, ys
end

function summarize_wins(xs, ys)
   breaks = sort(union(xs...))
   n = length(breaks)
   Z = zeros(length(xs), n)
   map(enumerate(zip(xs, ys))) do (k,(x, y))
       i = 1  # breaks index
       j = 1  # xs index
       while j <= length(x)
           # breaks is finer
           while i <= n && breaks[i] <= x[j]
               Z[k,i] = y[j]
               i += 1
           end
           j += 1
       end
   end
   return breaks, Z
end

theights(ts::TreeSequence) = theights(to_tskit(ts))

function theights(ts)
    xs = collect(ts.breakpoints())[1:end-1]
    th = map(ts.trees()) do tree
        length(tree.roots) > 1 ? NaN : tree.time(tree.root)
    end
    xs[2:end], th[1:end-1]
end

diffdiv(ts::TreeSequence, p1=1, p2=2; kwargs...) = diffdiv(to_tskit(ts), p1-1, p2-1; kwargs...)

function diffdiv(ts, pop1=0, pop2=1; windows=collect(ts.breakpoints()))
    ts.simplify(ts.samples())
    x0 = ts.samples(population=pop1)
    x1 = ts.samples(population=pop2)
    pi0 = ts.diversity(x0, mode="branch", windows=windows) ./ 2
    pi1 = ts.diversity(x1, mode="branch", windows=windows) ./ 2
    dxy = ts.divergence([x0, x1], mode="branch", windows=windows) ./ 2
    windows[2:end], pi0, pi1, dxy
end

function single_barrier_haploid(m, s, u, C, x, NA, NB)
    AA  = Architecture([BiAllelic(0.0)], [x])
    AB  = Architecture([BiAllelic(u  )], [x])
    MA  = GPMap([HaploidLocus(0.0, 1)])
    MB  = GPMap([HaploidLocus(-s , 1)])
    R   = LinearMap(C)
    nA = collect(1:NA)
    nB = collect(1:NB) .+ NA
    xA = [ ones(Int, 1) for _=1:NA]
    xB = [zeros(Int, 1) for _=1:NB]
    popA = WFPopulation(ploidy=Haploid(), N=NA, arch=AA, gpm=MA, recmap=R, x=xA, nodes=nA)
    popB = WFPopulation(ploidy=Haploid(), N=NB, arch=AB, gpm=MB, recmap=R, x=xB, nodes=nB)
    mpop = Fwd.TwoPopOneWay(m, popA, popB)
end

"""
    hmrecrate(xs::Vector)

Calculate harmonic mean recombination rate given a bunch of map positions.
"""
function hmrecrate(xs)
    rs  = Fwd.rec_matrix(xs)
    rhm = 0.0
    L = length(xs)
    for i=2:L
        for j=1:i-1
            rhm += 1/rs[i,j]
        end
    end
    (L*(L-1)/2)/rhm
end
