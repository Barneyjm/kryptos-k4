"""'Divided Berlin': K4 in two halves, each with its OWN periodic key (and optionally own alphabet).
Block 1 (pos 21-33) lies in the first half, block 2 (63-73) in the second, for any split point 34..63.
For periodic keys the phase is irrelevant, so each block is checked on its own."""
import itertools, numpy as np, k4
from round3 import ALPHA
B1 = list(range(21, 34)); B2 = list(range(63, 74)); CR = dict(zip(k4.CPOS.tolist(), k4.CPL))
def block_ok(block, PA, CA, sign, L):
    k = [(CA.index(k4.CT[i]) - sign * PA.index(CR[i])) % 26 for i in block]; d = {}
    for i, v in zip(block, k):
        if d.setdefault(i % L, v) != v: return -1
    return len(block) - len(d)              # number of checks the block passed
same, indep1, indep2 = [], [], []
for (pn, PA), (cn, CA) in itertools.product(ALPHA.items(), repeat=2):
    for sign in (1, -1):
        r1 = [(L, block_ok(B1, PA, CA, sign, L)) for L in range(1, 10)]; r1 = [(L, c) for L, c in r1 if c >= 4]
        r2 = [(L, block_ok(B2, PA, CA, sign, L)) for L in range(1, 8)]; r2 = [(L, c) for L, c in r2 if c >= 4]
        for L1, c1 in r1: indep1.append((c1, pn, cn, sign, L1))
        for L2, c2 in r2: indep2.append((c2, pn, cn, sign, L2))
        for (L1, c1), (L2, c2) in itertools.product(r1, r2): same.append((c1 + c2, pn, cn, sign, L1, L2))
print(f"alphabet pairs: {len(ALPHA)**2:,} x 2 cipher types")
print("East half (pos 22-34) fits a short key (>= 4 checks):", len(indep1), sorted(indep1, reverse=True)[:5])
print("West half (pos 64-74) fits a short key (>= 4 checks):", len(indep2), sorted(indep2, reverse=True)[:5])
print("BOTH halves, same alphabet + cipher type, own keys:", len(same), sorted(same, reverse=True)[:8])
# chance baseline: same procedure on shuffled ciphertext
rng = np.random.default_rng(0); real = k4.CT; cnt = []
for t in range(5):
    k4.CT = "".join(rng.permutation(list(real))); n = 0
    for (pn, PA), (cn, CA) in itertools.product(ALPHA.items(), repeat=2):
        for sign in (1, -1):
            a = any(block_ok(B1, PA, CA, sign, L) >= 4 for L in range(1, 10)); b = any(block_ok(B2, PA, CA, sign, L) >= 4 for L in range(1, 8))
            n += a and b
    cnt.append(n)
k4.CT = real; print("shuffled-ciphertext baseline, BOTH halves fitting:", cnt)
