
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

diffdiv(ts::TreeSequence, args...; kwargs...) = diffdiv(to_tskit(ts), args...; kwargs...)

function diffdiv(ts, pop1=0, pop2=1; windows=collect(ts.breakpoints()))
    ts.simplify(ts.samples())
    x0 = ts.samples(population=0)
    x1 = ts.samples(population=1)
    pi0 = ts.diversity(x0, mode="branch", windows=windows) ./ 2
    pi1 = ts.diversity(x1, mode="branch", windows=windows) ./ 2
    dxy = ts.divergence([x0, x1], mode="branch", windows=windows) ./ 2
    windows, pi0, pi1, dxy
end
