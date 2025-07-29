
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
