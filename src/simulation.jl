
simulate!(args...; kwargs...) = simulate!(default_rng(), args...; kwargs...)

function simulate!(rng::AbstractRNG, model, n; 
        callback=deepcopy, every=1, show_progress=true)
    pr = Progress(n; enabled=show_progress)
    nx = n ÷ every
    x0 = callback(model)
    xs = Vector{typeof(x0)}(undef, nx+1)
    j  = 1
    xs[j] = x0
    for i=1:n
        model = generation!(rng, model)
        next!(pr)
        if i % every == 0
            j += 1 
            xs[j] = callback(model)
        end
    end
    return xs
end

function winstat(stat, winsize, xs, ys)
    T = typeof(stat(ys))
    zs = []
    ywin = T[]
    win1 = winsize
    i = 1
    while i <= length(xs)
        while i <= length(xs) && xs[i] < win1 
            push!(ywin, ys[i])
            i += 1
        end
        win = ((win1-winsize), win1)
        push!(zs, (win, stat(ywin)))
        win1 += winsize
        ywin = T[]
    end
    return zs
end
