"""Linear complexity (Berlekamp-Massey over GF(2) and GF(13)) of the implied key runs (positions 21-33 and 63-73).
A key from ANY linear recurrence of order L (any coefficients/primer) mod 26 has linear complexity <= L mod 2 and mod 13."""
import itertools, numpy as np, random, k4
from round3 import ALPHA
def bm(s, p):
    s = [x % p for x in s]; C = [1]; B = [1]; L = 0; m = 1; b = 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1): d = (d + C[i] * s[n - i]) % p
        if d == 0: m += 1; continue
        coef = d * pow(b, -1, p) % p; T = C[:]
        C = C + [0] * (len(B) + m - len(C))
        for i, bi in enumerate(B): C[i + m] = (C[i + m] - coef * bi) % p
        if 2 * L <= n: L, B, b, m = n + 1 - L, T, d, 1
        else: m += 1
    return L, C[:L + 1]
def runs(k):
    k = dict(zip(k4.CPOS.tolist(), k)); return [k[i] for i in range(21, 34)], [k[i] for i in range(63, 74)]
def joint_ok(r1, r2, p, C):
    """does recurrence C (length L+1) generate both runs mod p?"""
    L = len(C) - 1
    for r in (r1, r2):
        for n in range(L, len(r)):
            if sum(C[i] * r[n - i] for i in range(L + 1)) % p: return False
    return True
def lc(k):
    r1, r2 = runs(k); out = []
    for p in (2, 13):
        La, Ca = bm(r1, p); Lb, Cb = bm(r2, p)
        out.append((max(La, Lb), joint_ok(r1, r2, p, Ca) or joint_ok(r1, r2, p, Cb)))
    return out
if __name__ == "__main__":
    CR = list(zip(k4.CPOS.tolist(), k4.CPL)); res = []
    for (pn, PA), (cn, CA) in itertools.product(ALPHA.items(), repeat=2):
        for kn, f in (("vig", lambda c, q: (c - q) % 26), ("beau", lambda c, q: (c + q) % 26)):
            k = [f(CA.index(k4.CT[i]), PA.index(p)) for i, p in CR]
            (l2, j2), (l13, j13) = lc(k)
            res.append((l13, l2, j13 and j2, pn, cn, kn))
    res.sort()
    print(f"alphabet pairs x types: {len(res):,}")
    print("lowest linear complexity found (mod 13, mod 2, one rule fits both runs):", res[:5])
    rnd = [lc([random.randrange(26) for _ in range(24)]) for _ in range(20000)]
    from collections import Counter
    print("random keys, mod-13 complexity distribution:", sorted(Counter(max(x[1][0], 0) for x in rnd).items()))
    print("K4 mod-13 complexity distribution:          ", sorted(Counter(r[0] for r in res).items()))
    # planted check: order-4 recurrence mod 26 with random coefficients
    random.seed(1); c = [random.randrange(26) for _ in range(4)]; seq = [random.randrange(26) for _ in range(4)]
    while len(seq) < 97: seq.append(sum(ci * seq[-1 - i] for i, ci in enumerate(c)) % 26)
    print("planted order-4 key ->", lc([seq[i] for i in k4.CPOS]))
