import itertools, numpy as np, random, k4
from linear import rank_consistent
from texts import KEYWORDS
CR = list(zip(k4.CPOS.tolist(), k4.CPL)); POS = np.array([i for i, _ in CR])
def unknown_split_test(ct, X, f):
    k = [f(X.index(ct[i]), X.index(p)) for i, p in CR]; out = []
    for L1 in range(2, 27):
        for L2 in range(L1 + 1, 27):
            rows = []
            for i, _ in CR:
                r = [0] * (L1 + L2); r[i % L1] = 1; r[L1 + i % L2] = 1; rows.append(r)
            a, r2 = rank_consistent(rows, k, 2); b, r13 = rank_consistent(rows, k, 13)
            if a and b and 24 - max(r2, r13) >= 4: out.append((24 - max(r2, r13), L1, L2))
    return out
# plant: PALIMPSEST + ABSCISSA (KA vigenere), random plaintext with cribs
random.seed(2); P = [random.choice(k4.AZ) for _ in range(97)]
for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
X = k4.KA; A, B = "PALIMPSEST", "ABSCISSA"
ct = "".join(X[(X.index(p) + X.index(A[(i + 3) % 10]) + X.index(B[(i + 5) % 8])) % 26] for i, p in enumerate(P))
print("planted (10, 8) found:", [x for x in unknown_split_test(ct, X, lambda c, p: (c - p) % 26) if x[1:] == (8, 10)])
# known clue words with all phase offsets, pairs, vectorized
words = sorted(set(KEYWORDS + ["SHADOW", "LUCID", "MEMORY", "POSITION", "QUAE", "DYAHR", "YAR", "WELTZEITUHR", "URANIA"]))
hits = []
for an, X in k4.ALPHS.items():
    cta = np.array([X.index(k4.CT[i]) for i in POS]); pta = np.array([X.index(p) for _, p in CR])
    for kn, need in (("vig", (cta - pta) % 26), ("beau", (cta + pta) % 26)):
        cache = {}
        for w in words:
            wa = np.array([X.index(c) for c in w]); L = len(w)
            cache[w] = np.array([wa[(POS + o) % L] for o in range(L)])     # (L offsets, 24)
        for w1, w2 in itertools.combinations(words, 2):
            for sgn in (1, -1):
                M = (cache[w1][:, None, :] + sgn * cache[w2][None, :, :]) % 26   # (o1, o2, 24)
                m = (M == need).sum(axis=2)
                for o1, o2 in zip(*np.nonzero(m >= 16)): hits.append((int(m[o1, o2]), an, kn, w1, w2, sgn, int(o1), int(o2)))
print("clue-word pairs with every starting offset, >= 16/24:", sorted(hits, reverse=True)[:8] or "none")
