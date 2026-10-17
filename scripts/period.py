import numpy as np, random, k4
random.seed(3)
def pic(a, L):
    tot = 0; pairs = 0
    for r in range(L):
        c = np.bincount(a[r::L], minlength=26); n = c.sum(); tot += (c*(c-1)).sum(); pairs += n*(n-1)
    return tot / pairs        # pooled coincidence rate within columns
A = np.array([ord(x)-65 for x in k4.CT])
sims = np.array([[pic(np.array([random.randrange(26) for _ in range(97)]), L) for L in range(1, 49)] for _ in range(20000)])
k = np.array([pic(A, L) for L in range(1, 49)])
print(" L   K4coinc  random mean   pct(K4 <= random)")
for L in range(1, 49):
    x = sims[:, L-1]; print(f"{L:2d}  {k[L-1]:.4f}   {x.mean():.4f}   {(x <= k[L-1]).mean():.3f}")
# pooled: sum of within-column coincidences across L=2..24 (correlated, use simulation)
stat = lambda row: row[1:24].mean()
print("mean coinc L=2..24: K4 %.4f, random pct %.4f" % (stat(k), (sims[:, 1:24].mean(axis=1) <= stat(k)).mean()))
