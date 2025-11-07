
n = 1000
ts0 = msprime.simulate(n, recombination_rate=1, random_seed=1)
ts1 = Fwd.from_tskit(ts0)

smpl = collect(1:10)
@code_warntype simplify(ts1, smpl)

@benchmark simplify(ts1, smpl)
# Time  (mean ± σ):   53.312 μs ± 663.323 μs  ┊ GC (mean ± σ):  16.91% ±  1.39%

# In [13]: %timeit ts0.simplify(range(10))
# 92 μs ± 127 ns per loop (mean ± std. dev. of 7 runs, 10,000 loops each)

# Other experiments (other n, larger sample), show that our implementation
# compares favorably against tskit (similar speed, sometimes faster, it
# seems). 
