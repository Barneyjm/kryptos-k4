import numpy as np, random, k4
random.seed(5)
def coinc_by_d(a):
    n = len(a); return np.array([(a[d:] == a[:-d]).sum() for d in range(1, n)])   # count of equal pairs at distance d
A = np.array([ord(x)-65 for x in k4.CT]); n = 97
c = coinc_by_d(A); pairs = np.array([n - d for d in range(1, n)])
print("total equal pairs:", c.sum(), " expected random:", round(pairs.sum()/26, 1))
sims = np.array([coinc_by_d(np.array([random.randrange(26) for _ in range(n)])) for _ in range(20000)])
ds = np.arange(1, n)
print(" m  residue: K4 equal pairs / expected-if-random   (pct of random sims <= K4)")
for m in (2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13):
    for r in range(m) if m <= 4 else (0,):
        sel = ds % m == r
        k = c[sel].sum(); e = pairs[sel].sum()/26; p = (sims[:, sel].sum(axis=1) <= k).mean()
        print(f"{m:2d}  d≡{r}: {k:3d} / {e:6.1f}   pct {p:.4f}")
print("distances of the equal pairs that are multiples of 3:", [int(d) for d in ds[(ds % 3 == 0)] for _ in range(c[d-1])][:40])
